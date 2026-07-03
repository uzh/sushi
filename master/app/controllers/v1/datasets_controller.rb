require 'json'

# Machine-callable, bearer-only registration API surface (v1 contract subset).
#
# Deliberately inherits from ActionController::Base (NOT ApplicationController):
# no Devise, no session, no CSRF token form. Authentication is fail-closed via a
# per-caller bearer token (ApiToken); authorization is least-privilege by the
# token's project scope. The manifest is received as content, never a path.
#
# Implemented: validate, register (idempotent), set_bfabric_id (set-once),
# destroy (deregister). Out of scope by design: audit hash-chain, attestation
# signing, reconciler, rate-limiting.
#
# Paths are /v1/... per the frozen v0.2 contract ([C4]); the registrar reaches
# this surface by base-URL swap (R-INV-10).
module V1
  class DatasetsController < ActionController::Base
    # No forgery protection cookie surface; this surface is bearer-only.
    protect_from_forgery with: :null_session

    # Single composable authorization gate (design v0.7). Ordered before_actions;
    # no action runs unless every step passes. Endpoint authority is checked
    # before the (network-bound) membership resolution, so a forbidden action
    # fails deterministically (403) rather than with a resolver-dependent 503.
    before_action :enforce_tls
    before_action :require_bearer_token
    before_action :reject_blank_user_login!
    before_action :authorize_endpoint!

    # A user token whose live membership resolver is unreachable fails closed with
    # 503 (a retryable availability failure, distinct from a 403 authz denial).
    rescue_from ApiToken::ResolverUnavailable do
      render json: { error: "authorization backend unavailable" }, status: :service_unavailable
    end

    # POST /v1/datasets/validate
    def validate
      return unless authorize_project!(body_param(:project_number))

      result = DatasetRegistrationService.validate(**manifest_args)
      render json: result, status: :ok
    end

    # POST /v1/datasets/register
    def register
      return unless authorize_project!(body_param(:project_number))

      result = DatasetRegistrationService.register(**manifest_args)
      render json: result[:body], status: result[:http]
    end

    # PUT /v1/datasets/:id/bfabric-id
    def set_bfabric_id
      data_set = find_scoped_data_set(params[:id])
      return unless data_set

      bfabric_id = body_param(:bfabric_id)
      if bfabric_id.to_s.empty?
        render json: { error: "bfabric_id is required" }, status: :unprocessable_entity
        return
      end

      result = DatasetRegistrationService.set_bfabric_id(data_set, bfabric_id)
      render json: result[:body], status: result[:http]
    end

    # DELETE /v1/datasets/:id
    def destroy
      data_set = DataSet.find_by(id: params[:id])
      # Idempotent: already-gone is a success.
      unless data_set
        render json: { ok: true, state: "COMPENSATED" }, status: :ok
        return
      end
      return unless scope_ok?(data_set)

      result = DatasetRegistrationService.deregister(data_set)
      render json: result[:body], status: result[:http]
    end

    private

    # --- filters ---------------------------------------------------------

    # INV-2T: confidential transport. Opt-in: enforced only when
    # SUSHI_API_ENFORCE_TLS=1 (set this where a TLS-terminating proxy fronts the
    # API). Off by default because the intended deployment reaches this surface
    # over an internal cleartext port (e.g. registrar -> http://fgcz-h-082:4071),
    # which a blanket production-https rule would wrongly reject.
    def enforce_tls
      return unless ENV["SUSHI_API_ENFORCE_TLS"] == "1"
      return if request.ssl? || request.headers["X-Forwarded-Proto"].to_s.downcase == "https"
      render json: { error: "TLS required" }, status: :forbidden
    end

    # INV-1 / INV-9: fail-closed bearer-only authentication.
    def require_bearer_token
      @api_token = ApiToken.authenticate(bearer_token)
      return if @api_token
      render json: { error: "unauthorized" }, status: :unauthorized
    end

    # A user token must carry a non-blank login; reject before any resolver call
    # (design v0.7 step 3). Static tokens are unaffected.
    def reject_blank_user_login!
      return unless @api_token&.user?
      return unless @api_token.login.to_s.strip.empty?
      render json: { error: "unauthorized" }, status: :unauthorized
    end

    # Deny-unless-listed endpoint authority for user tokens (design v0.7 step 4,
    # P-INV-4). A user principal may only reach the non-destructive whitelist;
    # anything else (today `destroy`/deregister; any future route) is 403. Static
    # principals pass unconditionally (P-INV-8).
    USER_ALLOWED_ACTIONS = %w[validate register set_bfabric_id].freeze
    def authorize_endpoint!
      return unless @api_token&.user?
      return if USER_ALLOWED_ACTIONS.include?(action_name)
      render json: { error: "action not permitted for this token" }, status: :forbidden
    end

    def bearer_token
      header = request.headers["Authorization"].to_s
      header[/\ABearer\s+(.+)\z/i, 1]
    end

    # --- authz helpers ---------------------------------------------------

    # Principal-aware membership test (design v0.7 step 6). Static → the token's
    # stored scope (unchanged). User → membership in the live-resolved set (W=0);
    # the set is resolved once per request and may raise ResolverUnavailable
    # (→ 503, handled by rescue_from).
    def token_allows_project?(project_number)
      return false if project_number.nil?
      if @api_token.user?
        (@user_allowed_projects ||= @api_token.allowed_projects).include?(project_number.to_i)
      else
        @api_token.in_scope?(project_number)
      end
    end

    # INV-8: request project_number must be authorized for the token.
    def authorize_project!(project_number)
      if project_number.to_s.empty?
        render json: { error: "project_number is required" }, status: :unprocessable_entity
        return false
      end
      unless token_allows_project?(project_number)
        render json: { error: "project #{project_number} out of scope" }, status: :forbidden
        return false
      end
      true
    end

    # IDOR guard: the target data_set's project must be authorized. Returns the
    # data_set, or nil after rendering. For a user token, a not-found dataset and
    # an out-of-scope existing dataset are BOTH 404 (indistinguishable) so dataset
    # IDs cannot be enumerated (design v0.7 step 6); static keeps 404/403 (P-INV-8).
    def find_scoped_data_set(id)
      data_set = DataSet.find_by(id: id)
      unless data_set
        render json: { error: "data_set #{id} not found" }, status: :not_found
        return nil
      end
      return nil unless scope_ok?(data_set)
      data_set
    end

    def scope_ok?(data_set)
      number = data_set.project&.number
      return true if token_allows_project?(number)

      if @api_token.user?
        render json: { error: "data_set #{data_set.id} not found" }, status: :not_found
      else
        render json: { error: "data_set out of scope" }, status: :forbidden
      end
      false
    end

    # --- request body ----------------------------------------------------

    def manifest_args
      {
        dataset_tsv:    body_param(:dataset_tsv),
        project_number: body_param(:project_number),
        name:           body_param(:name),
        parent_id:      body_param(:parent_id),
        order_id:       body_param(:order_id)
      }
    end

    def body_param(key)
      json_body[key.to_s]
    end

    def json_body
      @json_body ||= begin
        raw = request.body.read
        raw.to_s.empty? ? {} : JSON.parse(raw)
      rescue JSON::ParserError
        {}
      end
    end
  end
end

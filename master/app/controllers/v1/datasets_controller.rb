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

    before_action :enforce_tls
    before_action :require_bearer_token

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

    def bearer_token
      header = request.headers["Authorization"].to_s
      header[/\ABearer\s+(.+)\z/i, 1]
    end

    # --- authz helpers ---------------------------------------------------

    # INV-8: request project_number must be in the token's scope.
    def authorize_project!(project_number)
      if project_number.to_s.empty?
        render json: { error: "project_number is required" }, status: :unprocessable_entity
        return false
      end
      unless @api_token.in_scope?(project_number)
        render json: { error: "project #{project_number} out of scope" }, status: :forbidden
        return false
      end
      true
    end

    # IDOR guard: the target data_set's project must be in scope. Returns the
    # data_set, or nil after rendering 404/403.
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
      if number.nil? || !@api_token.in_scope?(number)
        render json: { error: "data_set out of scope" }, status: :forbidden
        return false
      end
      true
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

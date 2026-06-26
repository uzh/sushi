# Request specs for the bearer-only registration API
# (app/controllers/v1/datasets_controller.rb).
#
# NOTE: as of 2026-06-25 the legacy `test` environment does not boot under
# Rails 7 (config/environments/test.rb carries removed Rails-3 configs:
# `whiny_nils`, `mass_assignment_sanitizer`). Until that is repaired these
# specs cannot run; the behaviour they assert was verified end-to-end in the
# `development` environment via ActionDispatch::Integration::Session
# (auth matrix, idempotent replay, set-once bfabric-id, deregister, parent
# rejection). Keep this file so the suite is ready once the test env is fixed.

require 'spec_helper'

describe "Api::V1::Datasets", type: :request do
  let(:project_number) { 8888 }
  let(:tsv) { "Name\t[File]Read1\nsA\tp8888/data/a.fastq.gz\n" }
  let(:name) { "spec-#{SecureRandom.hex(4)}" }
  let!(:token_pair) { ApiToken.issue(name: "spec", scope: [project_number], ttl_days: 1) }
  let(:raw_token) { token_pair.first }
  let(:auth_headers) { { "Authorization" => "Bearer #{raw_token}", "CONTENT_TYPE" => "application/json" } }

  def body(overrides = {})
    { dataset_tsv: tsv, project_number: project_number, name: name }.merge(overrides).to_json
  end

  describe "authentication" do
    it "rejects a missing token with 401" do
      post "/v1/datasets/validate", params: body, headers: { "CONTENT_TYPE" => "application/json" }
      expect(response).to have_http_status(:unauthorized)
    end

    it "rejects an unknown token with 401" do
      post "/v1/datasets/validate", params: body,
           headers: { "Authorization" => "Bearer nope", "CONTENT_TYPE" => "application/json" }
      expect(response).to have_http_status(:unauthorized)
    end

    it "rejects an out-of-scope project with 403" do
      post "/v1/datasets/register", params: body(project_number: 7), headers: auth_headers
      expect(response).to have_http_status(:forbidden)
    end
  end

  describe "register" do
    it "registers and is idempotent on replay" do
      post "/v1/datasets/register", params: body, headers: auth_headers
      expect(response).to have_http_status(:ok)
      id = JSON.parse(response.body)["data_set_id"]
      expect(id).to be_present

      post "/v1/datasets/register", params: body, headers: auth_headers
      expect(JSON.parse(response.body)["data_set_id"]).to eq(id)
      expect(JSON.parse(response.body)["idempotent_replay"]).to be(true)
    end

    it "rejects a foreign/missing parent with 422" do
      post "/v1/datasets/register", params: body(parent_id: 99999999), headers: auth_headers
      expect(response).to have_http_status(:unprocessable_entity)
    end
  end

  describe "bfabric-id set-once" do
    it "accepts the same value and rejects a different one" do
      post "/v1/datasets/register", params: body, headers: auth_headers
      id = JSON.parse(response.body)["data_set_id"]

      put "/v1/datasets/#{id}/bfabric-id", params: { bfabric_id: 555 }.to_json, headers: auth_headers
      expect(response).to have_http_status(:ok)
      put "/v1/datasets/#{id}/bfabric-id", params: { bfabric_id: 555 }.to_json, headers: auth_headers
      expect(response).to have_http_status(:ok)
      put "/v1/datasets/#{id}/bfabric-id", params: { bfabric_id: 999 }.to_json, headers: auth_headers
      expect(response).to have_http_status(:conflict)
    end
  end

  describe "deregister" do
    it "is idempotent" do
      post "/v1/datasets/register", params: body, headers: auth_headers
      id = JSON.parse(response.body)["data_set_id"]

      delete "/v1/datasets/#{id}", headers: auth_headers
      expect(response).to have_http_status(:ok)
      delete "/v1/datasets/#{id}", headers: auth_headers
      expect(response).to have_http_status(:ok)
    end
  end
end

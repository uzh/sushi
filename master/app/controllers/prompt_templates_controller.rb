require 'net/http'
require 'json'

class PromptTemplatesController < ApplicationController
	before_action :authenticate_user!
	before_action :ensure_employee!

	PROMPT_TYPE = 'application_parameter_description'.freeze

	def show
		@template = fetch_template
		if @template.nil?
			flash.now[:alert] = 'Failed to load template from FastAPI.'
			@template = ''
		end
		render :application_parameter_description
	end

	def update
		new_template = params[:template].to_s
		if new_template.empty?
			flash[:alert] = 'Template cannot be empty.'
			@template = fetch_template
			return render :application_parameter_description
		end

		result = put_template(new_template)
		if result[:success]
			flash[:notice] = 'Template updated.'
			redirect_to '/prompt_templates/application_parameter_description'
		else
			flash.now[:alert] = "Failed to update template (status #{result[:status]} via #{result[:variant]}). #{result[:message]}"
			@template = new_template
			render :application_parameter_description
		end
	end

	def reset
		if reset_template
			flash[:notice] = 'Template reset to default.'
		else
			flash[:alert] = 'Failed to reset template.'
		end
		redirect_to '/prompt_templates/application_parameter_description'
	end

	private

	def ensure_employee!
		unless session[:employee]
			render plain: 'Forbidden', status: :forbidden
		end
	end

	def api_base
		ENV.fetch('PROMPT_API_BASE', 'http://fgcz-h-037:5002')
	end

	def fetch_template
		uri = URI.parse("#{api_base}/config/templates/#{PROMPT_TYPE}")
		req = Net::HTTP::Get.new(uri, { 'Accept' => 'application/json' })
		begin
			res = Net::HTTP.start(uri.hostname, uri.port, open_timeout: 2, read_timeout: 5) { |h| h.request(req) }
			return parse_template_response(res)
		rescue => e
			Rails.logger.warn("PromptTemplates#fetch_template failed: #{e.class}: #{e.message}")
			nil
		end
	end

	def parse_template_response(res)
		return nil unless res.is_a?(Net::HTTPSuccess)
		body = res.body.to_s
		begin
			json = JSON.parse(body)
			json['template'] || json['content'] || json['prompt'] || json['value'] || json['data'] || body
		rescue
			body
		end
	end

	def put_template(content)
		uri = URI.parse("#{api_base}/config/templates/#{PROMPT_TYPE}")
		attempts = []

		builders = [
			{
				variant: 'json.template',
				headers: { 'Content-Type' => 'application/json', 'Accept' => '*/*' },
				body: { template: content }.to_json
			},
			{
				variant: 'json.content',
				headers: { 'Content-Type' => 'application/json', 'Accept' => '*/*' },
				body: { content: content }.to_json
			},
			{
				variant: 'json.prompt',
				headers: { 'Content-Type' => 'application/json', 'Accept' => '*/*' },
				body: { prompt: content }.to_json
			},
			{
				variant: 'text.plain',
				headers: { 'Content-Type' => 'text/plain', 'Accept' => '*/*' },
				body: content
			},
			{
				variant: 'form.template',
				headers: { 'Content-Type' => 'application/x-www-form-urlencoded', 'Accept' => '*/*' },
				body: URI.encode_www_form(template: content)
			},
			{
				variant: 'form.content',
				headers: { 'Content-Type' => 'application/x-www-form-urlencoded', 'Accept' => '*/*' },
				body: URI.encode_www_form(content: content)
			}
		]

		builders.each do |b|
			begin
				req = Net::HTTP::Put.new(uri, b[:headers])
				req.body = b[:body]
				res = Net::HTTP.start(uri.hostname, uri.port, open_timeout: 2, read_timeout: 8) { |h| h.request(req) }
				code = res.code.to_s
				ok = res.is_a?(Net::HTTPSuccess)
				body_snippet = res.body.to_s[0, 160]
				attempts << { variant: b[:variant], status: code, ok: ok, body: body_snippet }
				if ok
					return { success: true, status: code, variant: b[:variant], message: '' }
				end
			rescue => e
				Rails.logger.warn("PromptTemplates#put_template attempt #{b[:variant]} failed: #{e.class}: #{e.message}")
				attempts << { variant: b[:variant], status: 'ERR', ok: false, body: e.message.to_s[0, 160] }
			end
		end

		last = attempts.last || { variant: 'n/a', status: 'n/a', body: '' }
		Rails.logger.warn("PromptTemplates#put_template all attempts failed: #{attempts.map { |a| "[#{a[:variant]} #{a[:status]}]" }.join(' ')}")
		{ success: false, status: last[:status], variant: last[:variant], message: last[:body] }
	end

	def reset_template
		uri = URI.parse("#{api_base}/config/templates/#{PROMPT_TYPE}/reset")
		req = Net::HTTP::Post.new(uri, { 'Accept' => 'application/json' })
		begin
			res = Net::HTTP.start(uri.hostname, uri.port, open_timeout: 2, read_timeout: 5) { |h| h.request(req) }
			res.is_a?(Net::HTTPSuccess)
		rescue => e
			Rails.logger.warn("PromptTemplates#reset_template failed: #{e.class}: #{e.message}")
			false
		end
	end
end



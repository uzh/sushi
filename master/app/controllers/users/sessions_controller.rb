# frozen_string_literal: true

class Users::SessionsController < Devise::SessionsController
  skip_before_action :verify_authenticity_token, only: [:create]
  # before_action :configure_sign_in_params, only: [:create]

  # GET /resource/sign_in
  def new
    super do
      flash[:alert] = nil if flash[:alert] == "You need to sign in or sign up before continuing."
    end
  end

  # POST /resource/sign_in
  # def create
  #   super
  # end

  # DELETE /resource/sign_out
  def destroy
    super do
      flash[:alert] = nil if flash[:alert] == "You need to sign in or sign up before continuing."
    end
  end

  # protected

  # If you have extra params to permit, append them to the sanitizer.
  # def configure_sign_in_params
  #   devise_parameter_sanitizer.permit(:sign_in, keys: [:attribute])
  # end
  def create
    self.resource = warden.authenticate!(auth_options)
    if resource
      set_flash_message!(:notice, :signed_in)
      sign_in(resource_name, resource)
      yield resource if block_given?
      respond_with resource, location: after_sign_in_path_for(resource)
    else
      set_flash_message!(:alert, :invalid)
      respond_to do |format|
        format.html { redirect_to new_user_session_path, alert: 'Failure of authentication' }
        format.js { render js: "alert('Failure of authentication');" }
      end
    end
  end

  protected

  def auth_options
    { scope: resource_name, recall: "#{controller_path}#new" }
  end
end

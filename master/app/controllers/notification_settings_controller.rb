class NotificationSettingsController < ApplicationController
  before_action :authenticate_user!
  before_action :check_permissions
  
  def show
    @notification_setting = current_user.notification_setting_or_create
  end
  
  def update
    @notification_setting = current_user.notification_setting_or_create
    
    if @notification_setting.update(notification_setting_params)
      flash[:notice] = 'Notification settings updated successfully.'
    else
      flash[:error] = 'Failed to update notification settings.'
    end
    
    redirect_to notification_setting_path
  end
  
  private
  
  def notification_setting_params
    params.require(:notification_setting).permit(:notification_enabled)
  end
  
  def check_permissions
    # Allow Employee users to manage their own notification settings
    unless session[:employee]
      flash[:error] = 'You do not have permission to manage notification settings.'
      redirect_to root_path
    end
  end
end 
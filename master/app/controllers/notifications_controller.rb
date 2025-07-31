class NotificationsController < ApplicationController
  before_action :authenticate_user!
  
  def index
    @notifications = current_user.recent_notifications(20)
  end
  
  def mark_as_read
    notification = current_user.notifications.find(params[:id])
    notification.mark_as_read!
    
    respond_to do |format|
      format.json { render json: { success: true } }
      format.html { redirect_back(fallback_location: notifications_path) }
    end
  end
  
  def mark_all_as_read
    current_user.notifications.unread.update_all(read: true)
    
    respond_to do |format|
      format.json { render json: { success: true } }
      format.html { redirect_back(fallback_location: notifications_path) }
    end
  end
  
  def unread_count
    count = current_user.unread_notifications_count
    render json: { count: count }
  end
end 
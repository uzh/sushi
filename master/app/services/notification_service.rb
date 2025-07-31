class NotificationService
  def initialize
    @logger = Rails.logger
  end
  
  # Send notification to employee users about dataset registration errors/warnings
  def notify_dataset_registration_issues(error_message = nil, warning_message = nil)
    employee_users = get_employee_users
    return if employee_users.empty?
    
    # Update notification settings for users who should receive notifications
    employee_users.each do |user|
      setting = user.notification_setting_or_create
      
      if error_message
        setting.update_error_date!
      elsif warning_message
        setting.update_warning_date!
      end
      
      # Only send notification if user has notifications enabled
      next unless setting.notification_enabled
      
      send_notification(user, error_message, warning_message)
      
      # Update last notification date after sending notification
      setting.update!(last_notification_date: Time.now)
    end
  end
  
  # Get all employee users with notifications enabled
  def get_employee_users
    User.joins(:notification_setting)
        .where(notification_settings: { notification_enabled: true })
        .select { |user| employee?(user) }
  end
  
  # Check if user is an employee
  def employee?(user)
    return false unless user
    
    # Use the same logic as in ApplicationHelper
    SushiFabric::Application.config.fgcz? && FGCZ.employee?(user.login)
  end
  
  # Send notification to a specific user
  def send_notification(user, error_message, warning_message)
    message = build_notification_message(error_message, warning_message)
    
    # Log the notification
    @logger.info "Dataset Registration Notification sent to #{user.login}: #{message}"
    
    # Create notification record in database
    begin
      notification_type = error_message ? 'dataset_registration_error' : 'dataset_registration_warning'
      Notification.create_dataset_registration_notification(user, message, notification_type)
      @logger.info "Notification record created for user #{user.login}"
    rescue => e
      @logger.error "Failed to create notification record for user #{user.login}: #{e.message}"
    end
    
    # Here you can implement additional notification delivery methods
    # For example: Send email notification
    # NotificationMailer.dataset_registration_alert(user, message).deliver_now
  end
  
  # Build notification message
  def build_notification_message(error_message, warning_message)
    timestamp = Time.now.strftime("%Y-%m-%d %H:%M:%S")
    
    if error_message
      "[ERROR] Dataset Registration Failed at #{timestamp}\n\nError: #{error_message}"
    elsif warning_message
      "[WARNING] Dataset Registration Warning at #{timestamp}\n\nWarning: #{warning_message}"
    else
      "[ERROR] Dataset Registration Failed at #{timestamp}\n\nAn error occurred during dataset registration process."
    end
  end
end 
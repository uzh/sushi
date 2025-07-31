class NotificationSetting < ActiveRecord::Base
  belongs_to :user
  
  validates :notification_enabled, inclusion: { in: [true, false] }
  
  # Check if there are new errors or warnings since last notification
  def has_new_errors_or_warnings?
    return false unless notification_enabled
    
    # Check if there are new errors since last notification
    if last_notification_date.nil? || last_error_date.nil?
      return false
    end
    
    last_error_date > last_notification_date || 
    (last_warning_date && last_warning_date > last_notification_date)
  end
  
  # Enable notifications and update last notification date
  def enable_notifications!
    update!(
      notification_enabled: true,
      last_notification_date: Time.now
    )
  end
  
  # Disable notifications
  def disable_notifications!
    update!(notification_enabled: false)
  end
  
  # Update error date and enable notifications if needed
  def update_error_date!
    update!(last_error_date: Time.now)
    enable_notifications! if !notification_enabled
  end
  
  # Update warning date and enable notifications if needed
  def update_warning_date!
    update!(last_warning_date: Time.now)
    enable_notifications! if !notification_enabled
  end
  
  # Update last notification date
  def update_last_notification_date!
    update!(last_notification_date: Time.now)
  end
  
  # Get or create notification setting for a user
  def self.for_user(user)
    find_or_create_by(user: user) do |setting|
      setting.notification_enabled = true
    end
  end
end

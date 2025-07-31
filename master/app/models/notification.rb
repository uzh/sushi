class Notification < ActiveRecord::Base
  belongs_to :user
  
  validates :message, presence: true
  validates :notification_type, presence: true
  
  scope :unread, -> { where(read: false) }
  scope :read, -> { where(read: true) }
  scope :recent, -> { order(created_at: :desc) }
  
  # Mark notification as read
  def mark_as_read!
    update!(read: true)
  end
  
  # Mark notification as unread
  def mark_as_unread!
    update!(read: false)
  end
  
  # Get unread notifications count for a user
  def self.unread_count_for_user(user)
    where(user: user, read: false).count
  end
  
  # Get recent notifications for a user
  def self.recent_for_user(user, limit = 10)
    where(user: user).recent.limit(limit)
  end
  
  # Create notification for dataset registration issues
  def self.create_dataset_registration_notification(user, message, notification_type = 'dataset_registration')
    create!(
      user: user,
      message: message,
      notification_type: notification_type,
      read: false
    )
  end
end

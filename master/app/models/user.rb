class User < ActiveRecord::Base
  # Include default devise modules. Others available are:
  # :token_authenticatable, :confirmable,
  # :lockable, :timeoutable and :omniauthable
  devise :ldap_authenticatable, :rememberable, :trackable

  # Setup accessible (or protected) attributes for your model
#  attr_accessible :login, :password, :password_confirmation, :remember_me, :selected_project
  # attr_accessible :title, :body
  has_many :data_sets
  has_one :notification_setting, dependent: :destroy
  has_many :notifications, dependent: :destroy
  
  # Get or create notification setting for this user
  def notification_setting_or_create
    notification_setting || create_notification_setting(notification_enabled: true)
  end
  
  # Get unread notifications count
  def unread_notifications_count
    Notification.unread_count_for_user(self)
  end
  
  # Get recent notifications
  def recent_notifications(limit = 10)
    Notification.recent_for_user(self, limit)
  end
end

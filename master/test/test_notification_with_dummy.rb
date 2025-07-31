#!/usr/bin/env ruby
# encoding: utf-8

# Test script for notification functionality with dummy error messages
ENV['RAILS_ENV'] = 'production'
ENV['DISABLE_DATABASE_ENVIRONMENT_CHECK'] = '1'

require 'bundler/setup'
require File.join(File.dirname(__FILE__), '..', 'config', 'environment')

puts "Testing notification functionality with dummy error messages..."

# Test notification service
notification_service = NotificationService.new

# Test 1: Dummy error notification
puts "\n1. Testing dummy error notification..."
error_message = "Test error: Dataset registration failed due to network timeout"
puts "  Sending error notification: #{error_message}"
begin
  notification_service.notify_dataset_registration_issues(error_message, nil)
  puts "  ✅ Error notification sent successfully"
rescue => e
  puts "  ❌ Error sending error notification: #{e.message}"
  puts "  Backtrace: #{e.backtrace.first(5).join("\n    ")}"
end

# Test 2: Dummy warning notification
puts "\n2. Testing dummy warning notification..."
warning_message = "Test warning: Some datasets were skipped due to missing metadata"
puts "  Sending warning notification: #{warning_message}"
begin
  notification_service.notify_dataset_registration_issues(nil, warning_message)
  puts "  ✅ Warning notification sent successfully"
rescue => e
  puts "  ❌ Error sending warning notification: #{e.message}"
  puts "  Backtrace: #{e.backtrace.first(5).join("\n    ")}"
end

# Test 3: Check notification records in database
puts "\n3. Checking notification records in database..."
begin
  total_notifications = Notification.count
  unread_notifications = Notification.unread.count
  
  puts "Total notifications: #{total_notifications}"
  puts "Unread notifications: #{unread_notifications}"
  
  # Show recent notifications
  puts "\nRecent notifications:"
  Notification.recent.limit(5).each do |notification|
    puts "  - #{notification.user.login}: #{notification.notification_type} (#{notification.read ? 'Read' : 'Unread'})"
    puts "    Message: #{notification.message[0..100]}..."
    puts "    Created: #{notification.created_at}"
    puts ""
  end
rescue => e
  puts "Database error: #{e.message}"
end

# Test 4: Check notification settings for first employee user
puts "\n4. Checking notification settings..."
begin
  # Find first employee user
  employee_users = User.all.select { |u| notification_service.employee?(u) }
  if employee_users.any?
    test_user = employee_users.first
    setting = test_user.notification_setting_or_create
    puts "User: #{test_user.login}"
    puts "Notification enabled: #{setting.notification_enabled}"
    puts "Last notification date: #{setting.last_notification_date}"
    puts "Last error date: #{setting.last_error_date}"
    puts "Last warning date: #{setting.last_warning_date}"
    puts "Unread notifications count: #{test_user.unread_notifications_count}"
    
    # Check if last_notification_date was updated
    if setting.last_notification_date
      puts "✅ Last notification date was updated successfully!"
    else
      puts "❌ Last notification date was not updated!"
    end
  else
    puts "No employee users found"
  end
rescue => e
  puts "Error checking notification settings: #{e.message}"
end

# Test 5: Test employee user detection
puts "\n5. Testing employee user detection..."
begin
  all_users = User.all
  employee_users = all_users.select { |u| notification_service.employee?(u) }
  puts "Total users: #{all_users.count}"
  puts "Employee users: #{employee_users.count}"
  employee_users.each do |user|
    puts "  - #{user.login}"
  end
rescue => e
  puts "Error detecting employee users: #{e.message}"
end

puts "\nDummy notification test completed!"
puts "\nNext steps:"
puts "1. Login as employee user to see notifications"
puts "2. Check the bell icon in the header"
puts "3. Click on notifications to mark as read"
puts "4. Visit /notifications to see all notifications"
puts "5. Check notification settings page for updated last notification date" 

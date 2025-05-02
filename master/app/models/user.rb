class User < ActiveRecord::Base
  # Include default devise modules. Others available are:
  # :token_authenticatable, :confirmable,
  # :lockable, :timeoutable and :omniauthable
  devise :ldap_authenticatable, :rememberable, :trackable

  # Setup accessible (or protected) attributes for your model
#  attr_accessible :login, :password, :password_confirmation, :remember_me, :selected_project
  # attr_accessible :title, :body
  has_many :data_sets

  def email
    "#{login}@example.com"  # fake to return email
  end
  def self.serialize_into_session(record)
    [record.id.to_s, nil]  # return two values (id, salt)
  end
  def self.serialize_from_session(*args)
    keys = args.first
    record = find_by(id: keys[0])
    record
  end
end

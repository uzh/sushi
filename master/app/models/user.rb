class User < ActiveRecord::Base
  # Include default devise modules. Others available are:
  # :token_authenticatable, :confirmable,
  # :lockable, :timeoutable and :omniauthable
  devise :rememberable, :trackable # :ldap_authenticatable temporarily disabled

  # Setup accessible (or protected) attributes for your model
#  attr_accessible :login, :password, :password_confirmation, :remember_me, :selected_project
  # attr_accessible :title, :body
  has_many :data_sets
end

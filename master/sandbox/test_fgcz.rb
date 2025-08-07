
require 'fgcz'
require 'net/ldap'

#load '/usr/local/ngseq/Repositories/gem-fgcz/lib/fgcz.rb'
puts FGCZ.get_user_groups 'tanguy'
puts "==="
puts FGCZ.get_user_projects 'tanguy'
puts '---'
puts FGCZ.get_user_groups('masaomi').class
print 'groups:'
p FGCZ.get_user_groups('masaomi')
puts FGCZ.get_user_projects('masaomi').class
print 'projects:'
p FGCZ.get_user_projects('masaomi')
puts '---'
p FGCZ.get_user_groups('hogehoge')
p FGCZ.get_user_projects('hogehoge')

puts "FGCZ.get_bioinformatician_users:"
p FGCZ.get_bioinformatician_users

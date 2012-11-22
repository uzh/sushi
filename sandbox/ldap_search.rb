#!/usr/local/ngseq/bin/ruby

# TODO:
# create a function authenticate(login, password) that returns an array of groups
# if the array is empty => cannot authenticate
# else use the groups to list available projects, samples, extracts...


require 'net/ldap'


user = "trxcopy"
password = "fXSC7o8g"


ldap_server = "fgcz-ldap.fgcz-net.unizh.ch"
ldap_base = "DC=FGCZ-NET,DC=unizh,DC=ch"
ldap_dn = "CN=#{user},OU=OU_Applications,OU=OU_Accounts,DC=FGCZ-NET,DC=unizh,DC=ch"
ldap_password = password


ldap = Net::LDAP.new
ldap.host = ldap_server
ldap.port = 389
ldap.auth ldap_dn, ldap_password
if not ldap.bind
    abort "Cannot bind, user not valid!"
end

filter = Net::LDAP::Filter.eq("cn", user)

ldap.search(:base => ldap_base, :filter => filter) do |entry|
    puts "#{user} is a member of:"
    entry.each do |attribute, values|
        if attribute == :memberof
            values.each do |value|
                group = value[/^CN=SG_([^,]+),/, 1]
                if group and group != ""
                    puts " - #{group}"
                end
            end
        end
    end
end

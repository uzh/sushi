#!/usr/local/ngseq/bin/ruby

require 'net/ldap'

ldap_server = "ldap://fgcz-ldap.fgcz-net.unizh.ch"
ldap_base = "DC=FGCZ-NET,DC=unizh,DC=ch"
ldap_dn = "CN=trxcopy,OU=OU_Applications,OU=OU_Accounts,DC=FGCZ-NET,DC=unizh,DC=ch"
ldap_password = "fXSC7o8g"


ldap = Net::LDAP.new
ldap.host = "fgcz-ldap.fgcz-net.unizh.ch"
ldap.port = 389
ldap.auth ldap_dn, ldap_password
if ldap.bind
    puts "OK"
else
    puts "ERR"
end

#!/usr/local/ngseq/bin/ruby

require 'net/ldap'

module FGCZ

    LDAP_SERVER = "fgcz-ldap.fgcz-net.unizh.ch"
    LDAP_BASE = "DC=FGCZ-NET,DC=unizh,DC=ch"
    LDAP_DN = "CN=${login},OU=OU_Applications,OU=OU_Accounts,DC=FGCZ-NET,DC=unizh,DC=ch"
    
=begin
    Returns the list of groups for a given user.
    If there is a problem while authenticating, the list is empty.
=end
    def FGCZ.get_user_groups(login, password)
        groups = Array.new

        my_dn = LDAP_DN
        my_dn["${login}"] = login

        ldap = Net::LDAP.new
        ldap.host = LDAP_SERVER
        ldap.port = 389
        ldap.auth my_dn, password
        if not ldap.bind
            puts "Cannot bind, user not valid!"
            return groups
        end

        filter = Net::LDAP::Filter.eq("cn", login)
        ldap.search(:base => LDAP_BASE, :filter => filter) do |entry|
            entry.each do |attribute, values|
                if attribute == :memberof
                    values.each do |value|
                        group = value[/^CN=SG_([^,]+),/, 1]
                        if group and group != ""
                            groups << group
                        end
                    end
                end
            end
        end

        return groups
    end
end

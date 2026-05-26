#!/usr/local/ngseq/bin/ruby

require 'savon'

client = Savon.client("http://fgcz-bfabric-demo.uzh.ch/bfabric/externaljob?wsdl")
puts client.wsdl.soap_actions

response = client.request :read do
    soap.body = {
        :parameters => {
            :login => ENV.fetch("BFABRIC_DEMO_LOGIN", "bemployee"),
            :password => ENV.fetch("BFABRIC_DEMO_PASSWORD"),
            :query => {
                :id => 3506
            }
        }
    }
end
puts response.to_hash


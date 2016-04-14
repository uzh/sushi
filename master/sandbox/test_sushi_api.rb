require 'json'
require 'net/http'

uri = URI('http://fgcz-s-034.uzh.ch:4000/api/index')
info = Hash.new
info[:project] = '1001'
info[:name] = 'DataSet Name'
info[:path] = '/srv/gstore/projects/p1001/HiSeq_20010101/dataset.tsv'

Net::HTTP.start(uri.host, uri.port) do |http|
    res = http.post uri.path, JSON.generate(info)
    puts res.code
    puts JSON.parse res.body
end




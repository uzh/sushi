import sys
import json
import urllib
import httplib
import socket

# Connect to the server
h = httplib.HTTPConnection("fgcz-s-034.uzh.ch", 4000, timeout=2)
try:
    h.connect()
except socket.error:
    print "Cannot contact the server!"
    sys.exit(1)

# POST data to the API
data = {'project': '1001', 'name': 'my name',
        'path': '/srv/gstore/projects/p1001/HiSeq_20010101/dataset.tsv'}
h.request("POST", "/api/index", json.dumps(data), {})

# Get the response
try:
    r = h.getresponse()
except socket.timeout:
    print "Connection timed out!"
    sys.exit(1)

headers = {}
for item in r.getheaders():
    headers[item[0].title()] = item[1]

# Print the response
print "Status:", r.status
#print "Headers:", headers
print "Body:", r.read()

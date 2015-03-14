import urllib
import urllib2
import ast
import json


# Read in the json for the mols
in_mols = ast.literal_eval(open("mols.json").read())
# Read in the json for the smis
in_smis = ast.literal_eval(open("smis.json").read())
print "Running smiles cluster test"
# Now set the url - this changes each time you run
url = 'http://127.0.0.1:8000/rdkit_cluster/cluster/'
# Now set the values in the get request - this is a simple JSON
values =  {'THRESHOLD' : '0.5', # The threshold to find similar molecules
          'FP_METHOD' : 'morgan', # The method to use
          'SIM_METHOD' : 'tanimoto', # The similarity method to use
          'SCREEN_LIB': in_smis }# The library to screen - now JSON molecule objects
# Make the json
data = json.dumps(values)
# Make the request object
req = urllib2.Request(url, data, {'Content-Type': 'application/json'})
# Now make the url
f = urllib2.urlopen(req)
# Read the response
response = f.read()
# Print the response
try:
    if max([x["values"]["cluster"]for x in ast.literal_eval(response)]) == 633:
        print "TEST PASSED"
    else:
        print max([x["values"]["cluster"]for x in ast.literal_eval(response)])
        print "TEST FAILED"
except:
    print "TEST FAILED"
# Now do the molstest
print "Running the mols cluster test"
# Now set the values in the get request - this is a simple JSON
values = {'THRESHOLD' : '0.5', # The threshold to find similar molecules
          'FP_METHOD' : 'morgan', # The method to use
          'SIM_METHOD' : 'tanimoto', # The similarity method to use
          'SCREEN_LIB': in_mols # The library to screen - now JSON molecule objects
 }
# Set the json
data = json.dumps(values)
# Make the request object
req = urllib2.Request(url, data, {'Content-Type': 'application/json'})
# Now make the url
f = urllib2.urlopen(req)
# Read the response
response = f.read()
# Print the response
try:
    if max([x["values"]["cluster"]for x in ast.literal_eval(response)]) == 633:
        print "TEST PASSED"
    else:
        print "TEST FAILED"
except:
    print "TEST FAILED"

print "Running the mols cluster test using a URL"
url = 'http://127.0.0.1:8000/rdkit_cluster/cluster_mol_body/'
# Now set the values in the get request - this is a simple JSON
values = {'THRESHOLD' : '0.5', # The threshold to find similar molecules
          'FP_METHOD' : 'morgan', # The method to use
          'SIM_METHOD' : 'tanimoto', # The similarity method to use
          'SCREEN_LIB': "http://demos.informaticsmatters.com/rest/files/json/31" # The library to screen - now JSON molecule objects
 }
# Set the put
data = urllib.urlencode(values)
# Make the request object
req = urllib2.Request(url, data)
# Now make the url
f = urllib2.urlopen(req)
# Read the response
response = f.read()
# Print the responsea
try:
    if max([x["values"]["cluster"]for x in ast.literal_eval(response)]) == 83:
        print "TEST PASSED"
    else:
        print "TEST FAILED"
except:
    print "TEST FAILED"



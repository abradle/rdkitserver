import urllib
import urllib2
import ast
import json
import sys
import os

my_plat = sys.platform
if my_plat == "darwin":
    docker_ip = str(os.environ["DOCKER_IP"])
else:
    docker_ip = "127.0.0.1"


print "Running mols screen test"
# Read in the json for the mols
in_mols = ast.literal_eval(open("mols.json").read())
# Read in the json for the smis
in_smis = ast.literal_eval(open("smis.json").read())
# Now set the url - this changes each time you run
url = 'http://'+docker_ip+':8000/rdkit_screen/screen/'
# Now set the values in the get request - this is a simple JSON
####### THIS IS WHAT YOU NEED TO SET
values = {'SMILES' : 'CCCCCC', #The smiles of your query mol. "."
          'THRESHOLD' : '0.5', # The threshold to find similar molecules
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
# Now run the test
try:
    if len(ast.literal_eval(response)["SCREEN 0"]["OUT_MOLS"]) == 4:
        print "TEST PASSED"
    else:
        print "TEST FAILED"
except:
    print "TEST FAILED"

print "Running smiles screen test"
# Now do the smis test
# Now set the values in the get request - this is a simple JSON
values = {'SMILES' : 'CCCCCC',#The smiles to screen againast
          'THRESHOLD' : '0.5', # The threshold to find similar molecules
          'FP_METHOD' : 'morgan', # The method to use
          'SIM_METHOD' : 'tanimoto', # The similarity method to use
          'SCREEN_LIB': in_smis # The library to screen - now JSON molecule objects
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
    if len(ast.literal_eval(response)["SCREEN 0"]["OUT_MOLS"]) == 4:
        print "TEST PASSED"
    else:
        print "TEST FAILED"
except:
    print "TEST FAILED"


print "Running json screen test"
url = 'http://'+docker_ip+':8000/rdkit_screen/screen_mol_body/'
# Now do the smis test
# Now set the values in the get request - this is a simple JSON
values = {"MOLS": "CCCCCCC",
          'THRESHOLD' : '0.2', # The threshold to find similar molecules
          'FP_METHOD' : 'morgan', # The method to use
          'SIM_METHOD' : 'tanimoto', # The similarity method to use
          'SCREEN_LIB': "http://demos.informaticsmatters.com/rest/files/json/31" 
}
# Set the json
# Set the put
data = urllib.urlencode(values)
# Make the request object
req = urllib2.Request(url, data)
# Now make the url
f = urllib2.urlopen(req)
# Read the response
response = f.read()
# Print the response
try:
    if len(ast.literal_eval(response)["SCREEN 0"]["OUT_MOLS"]) == 4:
        print "TEST PASSED"
    else:
        print "TEST FAILED"
except:
    print "TEST FAILED"

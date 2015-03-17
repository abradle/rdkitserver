sudo docker pull abradle/rdkitserver
OUTPUT="$(sudo docker run -d=true -p 8000:8000 -e "SERVER_NAME=localhost" -e "SERVER_PORT=8000" abradle/rdkitserver)"
echo STARTED SERVER
echo WAITING TO CONNECT
echo 1
sleep 1
echo 2
sleep 1
echo 3
sleep 1
echo 4
sleep 1
echo 5
sleep 1
curl -X POST --data-urlencode @smis.json "http://127.0.0.1:8000/rdkit_cluster/cluster_simple/?threshold=0.5&fp_method=morgan&sim_method=tanimoto" > out.test_one
cmp --silent out.test_one out.test_one_check && echo "TEST PASSED" || echo "TEST FAILED"
rm out.test_one
curl -X POST --data-urlencode @smis.json "http://127.0.0.1:8000/rdkit_screen/screen_simple/?smiles=CCCCCC&threshold=0.5&fp_method=morgan&sim_method=tanimoto" > out.test_two
cmp --silent out.test_two out.test_two_check && echo "TEST PASSED" || echo "TEST FAILED"
rm out.test_two
sudo docker stop $OUTPUT
sudo docker rm $OUTPUT

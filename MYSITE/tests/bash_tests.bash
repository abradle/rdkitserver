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
echo "RUNNING TEST ONE -> CLUSTERING"
curl -X POST --data-urlencode @smis.json "http://${DOCKER_IP:-"127.0.0.1"}:8000/rdkit_cluster/cluster_simple/?threshold=0.5&fp_method=morgan&sim_method=tanimoto" > out.test_one
cmp --silent out.test_one out.test_one_check && echo "TEST PASSED" || echo "TEST FAILED"
rm out.test_one
echo "RUNNING TEST TWO -> SCREENING"
curl -X POST --data-urlencode @smis.json "http://${DOCKER_IP:-"127.0.0.1"}:8000/rdkit_screen/screen_simple/?smiles=CCCCCC&threshold=0.5&fp_method=morgan&sim_method=tanimoto" > out.test_two
cmp --silent out.test_two out.test_two_check && echo "TEST PASSED" || echo "TEST FAILED"
rm out.test_two
echo "RUNNING TEST THREE -> CLUSTERING PREVIOUS ERROR"
curl -X POST --data-urlencode @db.json "http://${DOCKER_IP:-"127.0.0.1"}:8000/rdkit_cluster/cluster_simple/?threshold=0.5&fp_method=morgan&sim_method=tanimoto" > out.test_three
cmp --silent out.test_three out.test_three_check && echo "TEST PASSED" || echo "TEST FAILED"
rm out.test_three
sudo docker stop $OUTPUT
sudo docker rm $OUTPUT
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
echo "RUNNING TEST ONE -> MACCS and DICE"
curl -X POST -H "Content-Type: application/json" --data-urlencode @smis.json "http://${DOCKER_IP:-"127.0.0.1"}:8000/rdkit_cluster/cluster_simple/?threshold=0.5&fp_method=maccs&sim_method=dice" > out.test_one
cmp --silent out.test_one out.new_test_one_check && echo "TEST PASSED" || echo "TEST FAILED"
rm out.test_one
echo "RUNNING TEST TWO -> RDKIT_TOPO and TVERSKY "
curl -X POST -H "Content-Type: application/json" --data-urlencode @smis.json "http://${DOCKER_IP:-"127.0.0.1"}:8000/rdkit_screen/screen_simple/?smiles=CCCCCC&threshold=0.5&fp_method=rdkit_topo&sim_method=tversky" > out.test_two
cmp --silent out.test_two out.new_test_two_check && echo "TEST PASSED" || echo "TEST FAILED"
rm out.test_two
echo "RUNNING TEST THREE -> ATOM_PAIRS and COSINE"
curl -X POST -H "Content-Type: application/json" --data-urlencode @db.json "http://${DOCKER_IP:-"127.0.0.1"}:8000/rdkit_cluster/cluster_simple/?threshold=0.5&fp_method=atom_pairs&sim_method=cosine" > out.test_three
cmp --silent out.test_three out.new_test_three_check && echo "TEST PASSED" || echo "TEST FAILED"
rm out.test_three
sudo docker stop $OUTPUT
sudo docker rm $OUTPUT

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
python cluster_tests.py
python screen_tests.py
# Now stop it and kill
sudo docker stop $OUTPUT
sudo docker rm $OUTPUT

source /root/.bashrc
cd $RDBASE
cd External/INCHI-API
bash download-inchi.sh
cd $RDBASE
mkdir build
cd build
cmake -D RDK_BUILD_INCHI_SUPPORT=ON ..
make
make install

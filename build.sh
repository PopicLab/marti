mkdir -p build
cd build || exit 1
cmake -DMARTI_BUILD_TESTS=1 -DCMAKE_BUILD_TYPE=RELEASE ../marti
make

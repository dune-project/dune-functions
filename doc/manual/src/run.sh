

################################################################################
echo "Test with gcc"
mkdir -p build_gcc
cd build_gcc

g++ -std=c++11 -O3 -funroll-loops ../integration-test.cc -o integration-test

rm timings*
./integration-test | tee timings1
./integration-test | tee timings2
./integration-test | tee timings3
./integration-test | tee timings4

cd ..



################################################################################
echo "Test with gcc pgo"
mkdir -p build_gcc_pgo
cd build_gcc_pgo

g++ -std=c++11 -O3 -funroll-loops ../integration-test.cc -o integration-test -fprofile-generate
echo "Generating pgo profile"
./integration-test
g++ -std=c++11 -O3 -funroll-loops ../integration-test.cc -o integration-test -fprofile-use

rm timings*
./integration-test | tee timings1
./integration-test | tee timings2
./integration-test | tee timings3
./integration-test | tee timings4

cd ..



################################################################################
echo "Test with clang"
mkdir -p build_clang
cd build_clang

clang++ -std=c++11 -O3 -funroll-loops ../integration-test.cc -o integration-test

rm timings*
./integration-test | tee timings1
./integration-test | tee timings2
./integration-test | tee timings3
./integration-test | tee timings4

cd ..

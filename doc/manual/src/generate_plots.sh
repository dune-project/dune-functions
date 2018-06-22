

################################################################################
echo "Generating plots for gcc"
cd build_gcc

python ../plot.py

cp timings.pgf ../timings_gcc.pgf
cd ..



################################################################################
echo "Generating plots for gcc pgo"
cd build_gcc_pgo

python ../plot.py

cp timings.pgf ../timings_gcc_pgo.pgf
cd ..



################################################################################
echo "Generating plots for clang"
cd build_clang

python ../plot.py

cp timings.pgf ../timings_clang.pgf
cd ..

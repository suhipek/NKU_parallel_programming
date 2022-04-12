# mat_sub.sh
# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
pssh -h $PBS_NODEFILE mkdir -p /home/s2010056/3_simd 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/3_simd/gauss.cpp /home/s2010056/simd 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/3_simd/gauss.dat /home/s2010056/simd 1>&2
cd /home/s2010056/simd
for i in {1..32}
do 
    g++ -march=native -DN=$[64*i] ./gauss.cpp -o ./gauss_test
    ./gauss_test >> ./gauss_timing_$timestr.csv
done

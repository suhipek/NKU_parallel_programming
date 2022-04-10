# mat_sub.sh
# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
pssh -h $PBS_NODEFILE mkdir -p /home/s2010056/3_simd 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/3_simd/gauss /home/s2010056/3_simd 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/3_simd/gauss.dat /home/s2010056/3_simd 1>&2
/home/s2010056/3_simd/gauss

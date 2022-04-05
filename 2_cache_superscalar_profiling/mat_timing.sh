# mat_sub.sh
# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
pssh -h $PBS_NODEFILE mkdir -p /home/s2010056/2_cache_superscalar_profiling 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/2_cache_superscalar_profiling/matrix_product.cpp /home/s2010056/2_cache_superscalar_profiling 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/2_cache_superscalar_profiling/mar_vec.dat /home/s2010056/2_cache_superscalar_profiling 1>&2

for i in {1..150}
do 
    g++ -DN=$[10*i] /home/s2010056/2_cache_superscalar_profiling/matrix_product.cpp -o /home/s2010056/2_cache_superscalar_profiling/matrix_product
    /home/s2010056/2_cache_superscalar_profiling/matrix_product
done
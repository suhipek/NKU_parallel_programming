# mat_sub.sh
# !/bin/sh
pssh -h $PBS_NODEFILE mkdir -p /home/s2010056/2_cache_superscalar_profiling l>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/2_cache_superscalar_profiling/matrix_product /home/s2010056/2_cache_superscalar_profiling l>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/2_cache_superscalar_profiling/mar_vec.dat /home/s2010056/2_cache_superscalar_profiling l>&2
/s2010056/2_cache_superscalar_profiling/matrix_product
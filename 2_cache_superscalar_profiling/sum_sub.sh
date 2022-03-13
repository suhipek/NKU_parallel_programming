# mat_sub.sh
# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
pssh -h $PBS_NODEFILE mkdir -p /home/s2010056/2_cache_superscalar_profiling 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/2_cache_superscalar_profiling/cumulative /home/s2010056/2_cache_superscalar_profiling 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/2_cache_superscalar_profiling/mar_vec.dat /home/s2010056/2_cache_superscalar_profiling 1>&2
/home/s2010056/2_cache_superscalar_profiling/cumulative
perf record -e L1-dcache-load-misses,L1-dcache-loads,L1-dcache-store-misses,L1-dcache-stores,branch-loads,branch-load-misses,instructions,cycles,stalled-cycles-backend,stalled-cycles-frontend,branch-misses,cache-misses,cache-references -g -o /home/s2010056/2_cache_superscalar_profiling/perf_$timestr.data /home/s2010056/2_cache_superscalar_profiling/cumulative
pscp -h $PBS_NODEFILE /home/s2010056/2_cache_superscalar_profiling/perf_$timestr.data /home/s2010056/NKU_parallel_programming/2_cache_superscalar_profiling/ 1>&2
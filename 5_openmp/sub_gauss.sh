# !/bin/sh
timestr=$(date +%m_%d_%H_%M)


pssh -h $PBS_NODEFILE mkdir -p /home/s2010056/NKUPP 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/5_openmp/gauss.cpp /home/s2010056/NKUPP/5_openmp 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/5_openmp/gauss.dat /home/s2010056/NKUPP 1>&2

echo "行数,串行,无优化,静态,动态(8),动态(32),递减" >>/home/s2010056/NKUPP/5_openmp/gauss_timing_$timestr.csv
for i in {1..32}; do
    echo -n $((128 * i))"," >>/home/s2010056/NKUPP/5_openmp/gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DNUM_THREADS=1 -DSEPR=\",\" -DN=$((128 * i)) /home/s2010056/NKUPP/5_openmp/gauss.cpp -o /home/s2010056/NKUPP/5_openmp/gauss_test
    /home/s2010056/NKUPP/5_openmp/gauss_test >>/home/s2010056/NKUPP/5_openmp/gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DOMP_NO_OPT -DSEPR=\",\" -DN=$((128 * i)) /home/s2010056/NKUPP/5_openmp/gauss.cpp -o /home/s2010056/NKUPP/5_openmp/gauss_test
    /home/s2010056/NKUPP/5_openmp/gauss_test >>/home/s2010056/NKUPP/5_openmp/gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DOPT_CLAUSE=schedule\(static,N/NUM_THREADS\) -DSEPR=\",\" -DN=$((128 * i)) /home/s2010056/NKUPP/5_openmp/gauss.cpp -o /home/s2010056/NKUPP/5_openmp/gauss_test
    /home/s2010056/NKUPP/5_openmp/gauss_test >>/home/s2010056/NKUPP/5_openmp/gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DOPT_CLAUSE=schedule\(dynamic,N/NUM_THREADS/8\) -DSEPR=\",\" -DN=$((128 * i)) /home/s2010056/NKUPP/5_openmp/gauss.cpp -o /home/s2010056/NKUPP/5_openmp/gauss_test
    /home/s2010056/NKUPP/5_openmp/gauss_test >>/home/s2010056/NKUPP/5_openmp/gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DOPT_CLAUSE=schedule\(dynamic,N/NUM_THREADS/32\) -DSEPR=\",\" -DN=$((128 * i)) /home/s2010056/NKUPP/5_openmp/gauss.cpp -o /home/s2010056/NKUPP/5_openmp/gauss_test
    /home/s2010056/NKUPP/5_openmp/gauss_test >>/home/s2010056/NKUPP/5_openmp/gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DOPT_CLAUSE=schedule\(guided\) -DSEPR=\",\" -DN=$((128 * i)) /home/s2010056/NKUPP/5_openmp/gauss.cpp -o /home/s2010056/NKUPP/5_openmp/gauss_test
    /home/s2010056/NKUPP/5_openmp/gauss_test >>/home/s2010056/NKUPP/5_openmp/gauss_timing_$timestr.csv
    echo "" >>/home/s2010056/NKUPP/5_openmp/gauss_timing_$timestr.csv
done

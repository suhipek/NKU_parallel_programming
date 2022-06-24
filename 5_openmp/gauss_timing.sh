# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
lastfinish=$(date +%m_%d_%H_%M_%S)
num_th=("1" "4" "8" "12" "16" "20")
echo "行数,串行,无优化,静态,动态(8),动态(32),递减" >>./gauss_timing_$timestr.csv
for i in {1..32}; do
    echo -n $((128 * i))"," >>./gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DNUM_THREADS=1 -DSEPR=\",\" -DN=$((128 * i)) ./gauss.cpp -o ./gauss_test
    ./gauss_test >>./gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DOMP_NO_OPT -DSEPR=\",\" -DN=$((128 * i)) ./gauss.cpp -o ./gauss_test
    ./gauss_test >>./gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DOPT_CLAUSE=schedule\(static,N/NUM_THREADS\) -DSEPR=\",\" -DN=$((128 * i)) ./gauss.cpp -o ./gauss_test
    ./gauss_test >>./gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DOPT_CLAUSE=schedule\(dynamic,N/NUM_THREADS/8\) -DSEPR=\",\" -DN=$((128 * i)) ./gauss.cpp -o ./gauss_test
    ./gauss_test >>./gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DOPT_CLAUSE=schedule\(dynamic,N/NUM_THREADS/32\) -DSEPR=\",\" -DN=$((128 * i)) ./gauss.cpp -o ./gauss_test
    ./gauss_test >>./gauss_timing_$timestr.csv
    g++ -O2 -march=native -w -pthread -fopenmp -DOPT_CLAUSE=schedule\(guided\) -DSEPR=\",\" -DN=$((128 * i)) ./gauss.cpp -o ./gauss_test
    ./gauss_test >>./gauss_timing_$timestr.csv
    echo "" >>./gauss_timing_$timestr.csv
done

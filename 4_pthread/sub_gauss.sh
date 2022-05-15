# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
lastfinish=$(date +%m_%d_%H_%M_%S)
num_th=("1" "4" "8" "12" "16" "20")

pssh -h $PBS_NODEFILE mkdir -p /home/s2010056/4_pthread 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/4_pthread/gauss.cpp /home/s2010056/4_pthread 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/4_pthread/gauss.dat /home/s2010056/4_pthread 1>&2


for i in {1..32}; do
    echo -n $((128 * i))"," >>/home/s2010056/4_pthread/gauss_timing_$timestr.csv
    for j in {0..5}; do
        g++ -O0 -march=native -w -pthread -DNUM_THREADS=$((${num_th[j]})) -DN=$((128 * i)) /home/s2010056/4_pthread/gauss.cpp -o /home/s2010056/4_pthread/gauss_test
        echo "$(${num_th[i]}) threads"
        echo "now calculating: "$((128 * i))
        echo "time start: "$timestr
        echo "last finish: "$lastfinish
        echo "time now: "$(date +%m_%d_%H_%M_%S)
        /home/s2010056/4_pthread/gauss_test >>/home/s2010056/4_pthread/gauss_timing_$timestr.csv
        lastfinish=$(date +%m_%d_%H_%M_%S)
    done
    echo "">>/home/s2010056/4_pthread/gauss_timing_$timestr.csv
done

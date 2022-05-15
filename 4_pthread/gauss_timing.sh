# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
lastfinish=$(date +%m_%d_%H_%M_%S)
num_th=("1" "4" "8" "12" "16" "20")
for i in {1..32}; do
    echo -n $((128 * i))"," >>./gauss_timing_$timestr.csv
    for j in {0..5}; do
        g++ -O0 -march=native -w -pthread -DNUM_THREADS=$((${num_th[j]})) -DN=$((128 * i)) ./gauss.cpp -o ./gauss_test
        echo "$(${num_th[i]}) threads"
        echo "now calculating: "$((128 * i))
        echo "time start: "$timestr
        echo "last finish: "$lastfinish
        echo "time now: "$(date +%m_%d_%H_%M_%S)
        ./gauss_test >>./gauss_timing_$timestr.csv
        lastfinish=$(date +%m_%d_%H_%M_%S)
    done
    echo "">>./gauss_timing_$timestr.csv
done

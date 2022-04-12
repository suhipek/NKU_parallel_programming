# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
lastfinish=$(date +%m_%d_%H_%M_%S)
for i in {1..32}
do 
    g++ -march=native -DN=$[64*i] ./gauss.cpp -o ./gauss_test
    echo "now calculating: "$[64*i]
    echo "time start: "$timestr
    echo "last finish: "$lastfinish
    echo "time now: "$(date +%m_%d_%H_%M_%S)
    ./gauss_test >> ./gauss_timing_$timestr.csv
done
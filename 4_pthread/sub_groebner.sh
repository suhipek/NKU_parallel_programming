timestr=$(date +%m_%d_%H_%M)
pssh -h $PBS_NODEFILE mkdir -p /home/s2010056/4_pthread 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/4_pthread/groebner.cpp /home/s2010056/4_pthread 1>&2
data_path="/home/data/Groebner/"
num_th=("1" "4" "8" "12" "16" "20")

for file in $(ls ${data_path}); do
    if [ "$file" == "README.txt" ]; then
        continue
    else
        attr=(${file//_/ })
        # echo ${arrIN[1]}
        if [ "${attr[0]}" -gt "8" ]; then
            continue
        fi
        echo "${data_path}${file}/"
        for j in {0..5}; do
            echo ${num_th[j]}"threads"
            g++ -march=native -w -pthread -DDATA=\"${data_path}${file}/\" \
                -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} -DNUM_THREADS=$((${num_th[j]})) \
                /home/s2010056/4_pthread/groebner.cpp -o /home/s2010056/4_pthread/groebner
            /home/s2010056/4_pthread/groebner >>/home/s2010056/4_pthread/groebner_$timestr.csv
        done
        echo '' >>/home/s2010056/4_pthread/groebner_$timestr.csv
    fi
done

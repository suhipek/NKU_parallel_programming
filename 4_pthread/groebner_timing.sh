# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
data_path="../Groebner/"
num_th=("1" "4" "8" "12" "16" "20")

for file in $(ls ${data_path}); do
    if [ "$file" == "README.txt" ]; then
        continue
    else
        attr=(${file//_/ })
        # echo ${arrIN[1]}
        if [ "${attr[0]}" == "8" ]; then
            break
        fi
        echo "${data_path}${file}/"
        for j in {0..5}; do
            echo ${num_th[j]}"threads"
            g++ -march=native -w -pthread -DDATA=\"${data_path}${file}/\" \
                -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} \
                -DNUM_THREADS=$((${num_th[j]})) \
                ./groebner.cpp -o ./groebner
            ./groebner >>groebner_$timestr.csv
        done
        echo '' >>groebner_$timestr.csv
    fi
done

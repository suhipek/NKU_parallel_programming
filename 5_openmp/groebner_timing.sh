# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
data_path="../Groebner/"

echo ",串行,静态,动态(8),动态(32),递减" >>./gauss_timing_$timestr.csv

for file in $(ls ${data_path}); do
    if [ "$file" == "README.txt" ]; then
        continue
    else
        attr=(${file//_/ })
        # echo ${arrIN[1]}
        if [ "${attr[0]}" -lt "10" ]; then
            continue
        fi
        echo "${data_path}${file}/"
        echo -n "${file}," >>./groebner_$timestr.csv
        g++ -march=native -w -pthread -DDATA=\"${data_path}${file}/\" \
            -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} \
            -O2 -fopenmp -DSEPR=\",\" -DNUM_THREADS=1 \
            ./groebner.cpp -o ./groebner
        ./groebner >>groebner_$timestr.csv
        g++ -march=native -w -pthread -DDATA=\"${data_path}${file}/\" \
            -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} \
            -O2 -fopenmp -DSEPR=\",\" \
            -DOPT_CLAUSE=schedule\(static,ROW/NUM_THREADS\) \
            ./groebner.cpp -o ./groebner
        ./groebner >>groebner_$timestr.csv
        g++ -march=native -w -pthread -DDATA=\"${data_path}${file}/\" \
            -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} \
            -O2 -fopenmp -DSEPR=\",\" \
            -DOPT_CLAUSE=schedule\(dynamic,ROW/NUM_THREADS/8\) \
            ./groebner.cpp -o ./groebner
        ./groebner >>groebner_$timestr.csv
        g++ -march=native -w -pthread -DDATA=\"${data_path}${file}/\" \
            -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} \
            -O2 -fopenmp -DSEPR=\",\" \
            -DOPT_CLAUSE=schedule\(dynamic,ROW/NUM_THREADS/32\) \
            ./groebner.cpp -o ./groebner
        ./groebner >>groebner_$timestr.csv
        g++ -march=native -w -pthread -DDATA=\"${data_path}${file}/\" \
            -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} \
            -O2 -fopenmp -DSEPR=\",\" \
            -DOPT_CLAUSE=schedule\(guided\) \
            ./groebner.cpp -o ./groebner
        ./groebner >>groebner_$timestr.csv
        echo '' >>groebner_$timestr.csv
    fi
done

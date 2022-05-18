# mat_sub.sh
# !/bin/sh
timestr=$(date +%m_%d_%H_%M)
pssh -h $PBS_NODEFILE mkdir -p /home/s2010056/3_simd 1>&2
pscp -h $PBS_NODEFILE /home/s2010056/NKU_parallel_programming/3_simd/groebner.cpp /home/s2010056/3_simd 1>&2

data_path="/home/data/Groebner/"
cd /home/s2010056/3_simd
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
        g++ -march=native -w -DDATA=\"${data_path}${file}/\" \
            -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} \
            /home/s2010056/3_simd/groebner.cpp -o /home/s2010056/3_simd/groebner
        /home/s2010056/3_simd/groebner >>/home/s2010056/3_simd/groebner_$timestr.csv
        echo "unaligned finished"
        g++ -march=native -w -DALIGN -DDATA=\"${data_path}${file}/\" \
            -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} \
            /home/s2010056/3_simd/groebner.cpp -o /home/s2010056/3_simd/groebner
        /home/s2010056/3_simd/groebner >>/home/s2010056/3_simd/groebner_$timestr.csv
        echo '' >>/home/s2010056/3_simd/groebner_$timestr.csv
    fi
done
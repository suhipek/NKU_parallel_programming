# !/bin/sh
data_path="./Groebner/"
for file in $(ls ${data_path}); do
    if [ "$file" == "README.txt" ]; then
        continue
    else
        attr=(${file//_/ })
        # echo ${arrIN[1]}
        echo "${data_path}${file}/"
        
        g++ -march=native -w -DDATA=\"${data_path}${file}/\" \
            -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} \
            ./groebner.cpp -o ./groebner
        ./groebner >>groebner.csv
        echo "unaligned finished"
        g++ -march=native -w -DALIGN -DDATA=\"${data_path}${file}/\" \
            -DCOL=${attr[1]} -DELE=${attr[2]} -DROW=${attr[3]} \
            ./groebner.cpp -o ./groebner
        ./groebner >>groebner.csv
        echo '' >> groebner.csv
    fi
done

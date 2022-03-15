echo 'N,common_algo,optimized_algo'
for i in {1..90}
do 
    g++ -g -DN=$[16*i] ./matrix_product.cpp -o ./matrix_product
    ./matrix_product
done
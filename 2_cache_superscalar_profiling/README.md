# 第一次作业
先生成测试数据
```bash
g++ -O2 ./mar_vec_gen.cpp -o mar_vec_gen
./mar_vec_gen
```
- 矩阵列与向量的点积

    本质：cache优化

    ```bash
    g++ ./matrix_product.cpp -o matrix_product
    ./matrix_product
    ```
- 累加

    本质：乱序执行/分支预测的优化
    ```bash
    python3 ./unroll_gen.py # 生成unroll_sum.cpp
    g++ ./cumulative.cpp ./unroll_sum.cpp -o cumulative
    ./mar_vec_gen
    
    ```
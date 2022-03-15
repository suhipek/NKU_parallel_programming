#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include "unroll_sum.h"

//#define USE_FIXED_N
#define N 65536
// N MUST BE A POWER OF 2
#define REPT 1536

using namespace std;

void test(int (*func)(int *, int), const char *msg, int *arr, int len)
{
    timespec start, end;
    double time_used = 0;
    //cout << "result: " << func(arr, len) << "    ";
    clock_gettime(CLOCK_REALTIME, &start);
    //for (int i = 0; i < REPT*(int)pow(2,(20-(int)(log2(len)))); i++)
    for (int i = 0; i < REPT; i++)
        func(arr, len);
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    cout << time_used << ',';
}

int recursive_algo(int *arr, int len)
{
    for (int m = len; m > 1; m /= 2)
        for (int i = 0; i < m / 2; i++)
            arr[i] = arr[i * 2] + arr[i * 2 + 1];
    return arr[0];
}

int main()
{
    ifstream data("mar_vec.dat", ios::in | ios::binary);
    int arr[N];
    data.read((char *)arr, N * sizeof(int));
    data.close();

    cout << "N,common_algo,recursive_algo,unroll_algo_2,unroll_algo_8,unroll_algo_32,unroll_algo_128,unroll_algo_512,unroll_algo_2048" << endl;
    for(int i = 32; i<N; i*=2)
    {
        cout << i << ',';
        test(common_algo, "common_algo", arr, i);
        test(recursive_algo, "recursive_algo", arr, i);
        test(unroll_algo_2, "unroll_algo_2", arr, i);
        test(unroll_algo_8, "unroll_algo_8", arr, i);
        test(unroll_algo_32, "unroll_algo_32", arr, i);
        if(i>=128)test(unroll_algo_128, "unroll_algo_128", arr, i);
        if(i>=512)test(unroll_algo_512, "unroll_algo_512", arr, i);
        if(i>=2048)test(unroll_algo_2048, "unroll_algo_2048", arr, i);
        cout << endl;
    }
    cout << N << ',';
    test(common_algo, "common_algo", arr, N);
    test(recursive_algo, "recursive_algo", arr, N);
    test(unroll_algo_2, "unroll_algo_2", arr, N);
    // test(unroll_algo_4, "unroll_algo_4", arr, N);
    test(unroll_algo_8, "unroll_algo_8", arr, N);
    // test(unroll_algo_16, "unroll_algo_16", arr, N);
    test(unroll_algo_32, "unroll_algo_32", arr, N);
    // test(unroll_algo_64, "unroll_algo_64", arr, N);
    test(unroll_algo_128, "unroll_algo_128", arr, N);
    // test(unroll_algo_256, "unroll_algo_256", arr, N);
    test(unroll_algo_512, "unroll_algo_512", arr, N);
    // test(unroll_algo_1024, "unroll_algo_1024", arr, N);
    test(unroll_algo_2048, "unroll_algo_2048", arr, N);
    // test(unroll_algo_4096, "unroll_algo_4096", arr, N);
    // test(unroll_algo_8192, "unroll_algo_8192", arr, N);
    

    // for (int i = 256; i < N; i *= 2)
    // {
    //     cout << "data length: " << i << endl;
    //     test(common_algo, "common_algo", arr, i);
    //     test(unroll_algo_2, "unroll_algo_2", arr, i);
    //     test(unroll_algo_4, "unroll_algo_4", arr, i);
    //     test(unroll_algo_8, "unroll_algo_8", arr, i);
    //     test(unroll_algo_16, "unroll_algo_16", arr, i);
    //     test(unroll_algo_32, "unroll_algo_32", arr, i);
    //     test(unroll_algo_64, "unroll_algo_64", arr, i);
    //     test(unroll_algo_128, "unroll_algo_128", arr, i);
    //     test(unroll_algo_256, "unroll_algo_256", arr, i);
    //     if (i >= 512)
    //         test(unroll_algo_512, "unroll_algo_512", arr, i);
    //     if (i >= 1024)
    //         test(unroll_algo_1024, "unroll_algo_1024", arr, i);
    //     if (i >= 2048)
    //         test(unroll_algo_2048, "unroll_algo_2048", arr, i);
    //     if (i >= 4096)
    //         test(unroll_algo_4096, "unroll_algo_4096", arr, i);
    //     if (i >= 8192)
    //         test(unroll_algo_8192, "unroll_algo_8192", arr, i);
    //     test(recursive_algo, "recursive_algo", arr, i);
    // }

    return 0;
}
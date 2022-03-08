#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>

#define USE_FIXED_N
#define N 1048576
#define REPT 500

using namespace std;

int common_algo(int *arr, int len)
{
    int total = 0;
    for (int i = 0; i < len; i++)
        total += arr[i];
    return total;
}

int dual_algo(int *arr, int len)
{
    int total1 = 0, total2 = 0;
    for (int i = 0; i < len; i += 2)
    {
        total1 += arr[i];
        total2 += arr[i + 1];
    }
    return total1 + total2;
}

int quad_algo(int *arr, int len)
{
    int total1 = 0, total2 = 0, total3 = 0, total4 = 0;
    for (int i = 0; i < len; i += 4)
    {
        total1 += arr[i];
        total2 += arr[i + 1];
        total3 += arr[i + 2];
        total4 += arr[i + 3];
    }
    return total1 + total2 + total3 + total4;
}

#define SUM             \
    total[m] += arr[j]; \
    j++;
#define SUM2 SUM SUM
#define SUM4 SUM2 SUM2
#define SUM8 SUM4 SUM4
#define SUM16 SUM8 SUM8
#define SUM32 SUM16 SUM16
#define SUM64 SUM32 SUM32
#define SUM128 SUM64 SUM64
#define SUM256 SUM128 SUM128
#define SUM512 SUM256 SUM256
#define SUM1024 SUM512 SUM512
#define SUM2048 SUM1024 SUM1024
#define SUM4096 SUM2048 SUM2048
#define SUM8192 SUM4096 SUM4096

int unroll_algo(int *arr, int len)
{
    int j, m, total[len / 8192] = {0};
    for (int i = 0; i < len - 8192; i += 8192)
    {
        j = i;
        m = i / 8192;
        SUM8192
    }
    return common_algo(total, len / 8192);
}

int main()
{
    ifstream data("mar_vec.dat", ios::in | ios::binary);
    int arr[N];
    data.read((char *)arr, N * sizeof(int));
    data.close();

    timespec start, end;
    double time_used = 0;
    clock_gettime(CLOCK_REALTIME, &start);
    for (int i = 0; i < REPT; i++)
        common_algo(arr, N);
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    cout << "common_algo: " << time_used << endl;

    time_used = 0;
    clock_gettime(CLOCK_REALTIME, &start);
    for (int i = 0; i < REPT; i++)
        dual_algo(arr, N);
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    cout << "dual_algo: " << time_used << endl;

    time_used = 0;
    clock_gettime(CLOCK_REALTIME, &start);
    for (int i = 0; i < REPT; i++)
        quad_algo(arr, N);
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    cout << "quad_algo: " << time_used << endl;

    time_used = 0;
    clock_gettime(CLOCK_REALTIME, &start);
    for (int i = 0; i < REPT; i++)
        unroll_algo(arr, N);
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    cout << "unroll_algo: " << time_used << endl;

    int ret_common = common_algo(arr, N);
    int ret_dual = dual_algo(arr, N);

    if (ret_dual != ret_common)
        cout << "results mismatch" << endl;
    return 0;
}
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>

#define USE_FIXED_N
#define N 1024
#define REPT 100

using namespace std;

int *common_algo(int mat[N][N], int vec[N], int n)
{

    int *sum = new int[n];
    for (int i = 0; i < n; i++)
    {
        sum[i] = 0;
        for (int j = 0; j < n; j++)
            sum[i] += mat[j][i] * vec[j];
    }
    return sum;
}

int *optimized_algo(int mat[N][N], int vec[N], int n)
{
    int *sum = new int[n];
    for (int i = 0; i < n; i++)
        sum[i] = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            sum[j] += mat[i][j] * vec[i];
    return sum;
}

int main()
{
    ifstream data("mar_vec.dat", ios::in | ios::binary);
    int mat[N][N];
    int vec[N];
    data.read((char *)mat, N * N * sizeof(int));
    data.read((char *)vec, N * sizeof(int));
    data.close();

    timespec start, end;
    double time_used = 0;
    clock_gettime(CLOCK_REALTIME, &start);
    for (int i = 0; i < REPT; i++)
    {
        int *ret_common = common_algo(mat, vec, N);
        delete[] ret_common;
    }
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    cout << "common_algo: " << time_used << endl;

    time_used = 0;
    clock_gettime(CLOCK_REALTIME, &start);
    for (int i = 0; i < REPT; i++)
    {
        int *ret_optimized = optimized_algo(mat, vec, N);
        delete[] ret_optimized;
    }
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    cout << "optimized_algo: " << time_used << endl;

    // for(int i=0;i<N;i++)
    //     if(ret_common[i] != ret_optimized[i])
    //     {
    //         cout << "results not matched!" << endl;
    //         break;
    //     }

    // delete[] ret_common;
    // delete[] ret_optimized;
    return 0;
}
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>

#define USE_FIXED_N
//#define N 1024
#define REPT 200

using namespace std;

int mat[N][N];
int vec[N];

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

void test(int *(*func)(int mat[N][N], int vec[N], int n),
          const char *msg, int mat[N][N], int vec[N], int n)
{
    timespec start, end;
    double time_used = 0;
    int *ret = func(mat, vec, n);
    // cout << "result[0]: " << *ret << "    ";
    delete[] ret;
    clock_gettime(CLOCK_REALTIME, &start);
    for (int i = 0; i < REPT; i++)
    {
        int *ret = func(mat, vec, N);
        delete[] ret;
    }
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    // cout << msg << ": " << time_used << endl;
    cout << time_used << ',';
}

int main()
{
    ifstream data("mar_vec.dat", ios::in | ios::binary);
    data.read((char *)mat, N * N * sizeof(int));
    data.read((char *)vec, N * sizeof(int));
    data.close();

    cout << N << ',';
    test(common_algo, "common_algo", mat, vec, N);
    test(optimized_algo, "optimized_algo", mat, vec, N);
    cout << endl;

    return 0;
}
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <cmath>

//#include <arm_neon.h>
#include "NEON_2_SSE.h"

#define REPT 100
#define N 256
#define ele_t float
//#define DEBUG

using namespace std;

ele_t mat[N][N];
ele_t new_mat[N][N];

void test(void (*func)(ele_t[N][N], int), const char *msg, ele_t mat[N][N], int len)
{
    timespec start, end;
    double time_used = 0;
    // cout << "result: " << func(arr, len) << "    ";
    clock_gettime(CLOCK_REALTIME, &start);
    // for (int i = 0; i < REPT*(int)pow(2,(20-(int)(log2(len)))); i++)
    for (int i = 0; i < REPT; i++)
        func(mat, len);
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    cout << time_used << ',';
}

void LU(ele_t mat[N][N], int n)
{
    //ele_t new_mat[N][N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            new_mat[i][j] = mat[i][j];

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            if (new_mat[i][i] == 0)
                continue;
            ele_t div = new_mat[j][i] / new_mat[i][i];
            // cout << new_mat[j][i] << '/' << new_mat[i][i] << '=' << div << endl;
            for (int k = i; k < n; k++)
                new_mat[j][k] -= new_mat[i][k] * div;
        }

#ifdef DEBUG
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << new_mat[i][j] << ' ';
        cout << endl;
    }
#endif
}

void LU_simd(ele_t mat[N][N], int n)
{
    //ele_t new_mat[N][N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            new_mat[i][j] = mat[i][j];

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            if (new_mat[i][i] == 0)
                continue;
            ele_t div = new_mat[j][i] / new_mat[i][i];
            float32x4_t div4 = vmovq_n_f32(-1 * div);
            // cout << new_mat[j][i] << '/' << new_mat[i][i] << '=' << div << endl;
            for (int k = i; k < n; k += 4)
            {
                float32x4_t mat_j = vld1q_f32(new_mat[j] + k);
                float32x4_t mat_i = vld1q_f32(new_mat[i] + k);
                float32x4_t res = vmlaq_f32(mat_j, div4, mat_i);
                vst1q_f32(new_mat[j] + k, res);
            }
        }

#ifdef DEBUG
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << new_mat[i][j] << ' ';
        cout << endl;
    }
#endif
}

int main()
{

    ifstream data("gauss.dat", ios::in | ios::binary);
    data.read((char *)mat, N * N * sizeof(ele_t));
    data.close();

#ifndef DEBUG
    test(LU, "commone algo:", mat, N);
    test(LU_simd, "optimized algo:", mat, N);
#else
    LU(mat, N);
    LU_simd(mat, N);
#endif

    return 0;
}
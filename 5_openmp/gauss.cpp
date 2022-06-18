#include <pthread.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <string.h>
#include <omp.h>

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

#ifdef __amd64__
#include <immintrin.h>
#include "../NEON_2_SSE.h"
#define USE_SSE4
#endif

#define ZERO (float)1e-5
#ifndef N
#define N 4096
#define NUM_THREADS 16
#define NUM_BLOCKS 32
#endif
#define REPT 1
#define ele_t float
//#define DEBUG

#ifndef SEPR
#define SEPR endl
#endif

using namespace std;

ele_t new_mat[N][N] __attribute__((aligned(64)));
ele_t mat[N][N];

void test(void (*func)(ele_t[N][N], int), const char *msg, ele_t mat[N][N], int len)
{
    timespec start, end;
    double time_used = 0;
    clock_gettime(CLOCK_REALTIME, &start);
    // for (int i = 0; i < REPT*(int)pow(2,(20-(int)(log2(len)))); i++)
    for (int i = 0; i < REPT; i++)
        func(mat, len);
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    cout << time_used << SEPR;
}

void LU(ele_t mat[N][N], int n)
{
    memcpy(new_mat, mat, sizeof(ele_t) * N * N);

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
    cout << endl;
#endif
}

void LU_simd(ele_t mat[N][N], int n)
{
    // ele_t new_mat[N][N];
    memcpy(new_mat, mat, sizeof(ele_t) * N * N);

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            if (new_mat[i][i] == 0)
                continue;
            ele_t div = new_mat[j][i] / new_mat[i][i];
            float32x4_t div4 = vmovq_n_f32(div);
            float32x4_t mat_j;
            float32x4_t mat_i;
            float32x4_t res;
            // cout << new_mat[j][i] << '/' << new_mat[i][i] << '=' << div << endl;
            for (int k = i / 4 * 4; k < n; k += 4)
            {
                mat_j = vld1q_f32(new_mat[j] + k);
                mat_i = vld1q_f32(new_mat[i] + k);
                res = vmlsq_f32(mat_j, div4, mat_i);
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
    cout << endl;
#endif
}

void LU_omp(ele_t mat[N][N], int n)
{
    // ele_t new_mat[N][N];
    memcpy(new_mat, mat, sizeof(ele_t) * N * N);
    int i = 0;

    for (i = 0; i < n; i++)
#pragma omp parallel for
        for (int j = i + 1; j < n; j++)
        {
            if (fabs(new_mat[i][i]) < ZERO)
                continue;
            ele_t div = new_mat[j][i] / new_mat[i][i];
            float32x4_t div4 = vmovq_n_f32(div);
            float32x4_t mat_j;
            float32x4_t mat_i;
            float32x4_t res;
            // cout << new_mat[j][i] << '/' << new_mat[i][i] << '=' << div << endl;
            for (int k = i / 4 * 4; k < n; k += 4)
            {
                mat_j = vld1q_f32(new_mat[j] + k);
                mat_i = vld1q_f32(new_mat[i] + k);
                res = vmlsq_f32(mat_j, div4, mat_i);
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
    cout << endl;
#endif
}

void LU_omp_opt(ele_t mat[N][N], int n)
{
    memcpy(new_mat, mat, sizeof(ele_t) * N * N);

#pragma omp parallel num_threads(NUM_THREADS)
    for (int i = 0; i < n; i++)
    {
#pragma omp for schedule(dynamic, N/NUM_THREADS/NUM_BLOCKS)
        for (int j = i + 1; j < n; j++)
        {
            if (fabs(new_mat[i][i]) < ZERO)
                continue;
            ele_t div = new_mat[j][i] / new_mat[i][i];
#pragma omp simd simdlen(4)
            for (int k = i; k < n; k++)
                new_mat[j][k] -= new_mat[i][k] * div;
        }
    }

#ifdef DEBUG
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << new_mat[i][j] << ' ';
        cout << endl;
    }
    cout << endl;
#endif
}

int main()
{
    ifstream data("../gauss.dat", ios::in | ios::binary);
    data.read((char *)mat, N * N * sizeof(ele_t));
    data.close();

    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //         cout << mat[i][j] << ' ';
    //     cout << endl;
    // }

#ifndef DEBUG
    // test(LU, "commone algo: ", mat, N);
    // test(LU, "NEON/SSE: ", mat, N);
    // test(LU_simd, "NEON/SSE: ", mat, N);
    test(LU_omp, "NEON/SSE: ", mat, N);
    test(LU_omp_opt, "NEON/SSE: ", mat, N);
#else
    cout << endl;
    LU(mat, N);
    cout << endl
         << endl;
    LU_omp_opt(mat, N);
    cout << endl
         << endl;
#endif
    return 0;
}

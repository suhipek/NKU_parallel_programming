#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <cmath>

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

#ifdef __amd64__
#include <immintrin.h>
#include "NEON_2_SSE.h"
#define USE_SSE4
#endif

#ifndef N
#define N 1024
#endif
#define REPT 1
#define ele_t float
//#define DEBUG

using namespace std;

ele_t new_mat[N][N] __attribute__((aligned(64)));
ele_t mat[N][N];

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
    // ele_t new_mat[N][N];
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
    cout << endl;
#endif
}

void LU_simd(ele_t mat[N][N], int n)
{
    // ele_t new_mat[N][N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            new_mat[i][j] = mat[i][j];

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
            for (int k = i; k < n; k += 4)
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

void LU_simd_Aligned(ele_t mat[N][N], int n)
{
    // ele_t new_mat[N][N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            new_mat[i][j] = mat[i][j];

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

#ifdef __amd64__
void LU_sse_fma(ele_t mat[N][N], int n)
{
    // ele_t new_mat[N][N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            new_mat[i][j] = mat[i][j];

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
                res = _mm_fnmadd_ps(mat_i, div4, mat_j);
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

void LU_avx(ele_t mat[N][N], int n)
{
    // ele_t new_mat[N][N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            new_mat[i][j] = mat[i][j];

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            if (new_mat[i][i] == 0)
                continue;
            ele_t div = new_mat[j][i] / new_mat[i][i];
            __m256 div8 = _mm256_set1_ps(div);
            __m256 mat_j;
            __m256 mat_i;
            __m256 res;
            // cout << new_mat[j][i] << '/' << new_mat[i][i] << '=' << div << endl;
            for (int k = i; k < n; k += 8)
            {
                mat_j = _mm256_loadu_ps(new_mat[j] + k);
                mat_i = _mm256_loadu_ps(new_mat[i] + k);
                res = _mm256_fnmadd_ps(mat_i, div8, mat_j);
                _mm256_storeu_ps(new_mat[j] + k, res);
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

void LU_avx_aligned(ele_t mat[N][N], int n)
{
    // ele_t new_mat[N][N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            new_mat[i][j] = mat[i][j];

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            if (new_mat[i][i] == 0)
                continue;
            ele_t div = new_mat[j][i] / new_mat[i][i];
            __m256 div8 = _mm256_set1_ps(div);
            __m256 mat_j;
            __m256 mat_i;
            __m256 res;
            // cout << new_mat[j][i] << '/' << new_mat[i][i] << '=' << div << endl;
            for (int k = i / 8 * 8; k < n; k += 8)
            {
                mat_j = _mm256_load_ps(new_mat[j] + k);
                mat_i = _mm256_load_ps(new_mat[i] + k);
                res = _mm256_fnmadd_ps(mat_i, div8, mat_j);
                _mm256_store_ps(new_mat[j] + k, res);
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
#endif

#ifdef __AVX512F__
void LU_avx512(ele_t mat[N][N], int n)
{
    // ele_t new_mat[N][N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            new_mat[i][j] = mat[i][j];

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            if (new_mat[i][i] == 0)
                continue;
            ele_t div = new_mat[j][i] / new_mat[i][i];
            __m512 div16 = _mm512_set1_ps(div);
            __m512 mat_j;
            __m512 mat_i;
            __m512 res;
            // cout << new_mat[j][i] << '/' << new_mat[i][i] << '=' << div << endl;
            for (int k = i; k < n; k += 16)
            {
                mat_j = _mm512_loadu_ps(new_mat[j] + k);
                mat_i = _mm512_loadu_ps(new_mat[i] + k);
                res = _mm512_fnmadd_ps(mat_i, div16, mat_j);
                _mm512_storeu_ps(new_mat[j] + k, res);
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

void LU_avx512_aligned(ele_t mat[N][N], int n)
{
    // ele_t new_mat[N][N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            new_mat[i][j] = mat[i][j];

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            if (new_mat[i][i] == 0)
                continue;
            ele_t div = new_mat[j][i] / new_mat[i][i];
            __m512 div16 = _mm512_set1_ps(div);
            __m512 mat_j;
            __m512 mat_i;
            __m512 res;
            // cout << new_mat[j][i] << '/' << new_mat[i][i] << '=' << div << endl;
            for (int k = i / 16 * 16; k < n; k += 16)
            {
                mat_j = _mm512_load_ps(new_mat[j] + k);
                mat_i = _mm512_load_ps(new_mat[i] + k);
                res = _mm512_fnmadd_ps(mat_i, div16, mat_j);
                _mm512_store_ps(new_mat[j] + k, res);
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
#endif

int main()
{

    ifstream data("gauss.dat", ios::in | ios::binary);
    data.read((char *)mat, N * N * sizeof(ele_t));
    data.close();
    cout << N << ',';

#ifdef NO_ALIGN_INFO
    cout << "new_mat addr:" << &new_mat << "    ";
    cout << "% 256 == " << (long long)new_mat % 256 << endl;
    cout << "alignof(new_mat): " << alignof(new_mat) << endl;
#endif

#ifndef DEBUG
    test(LU, "commone algo: ", mat, N);
    test(LU_simd, "NEON/SSE: ", mat, N);
    test(LU_simd_Aligned, "NEON/SSE Aligned: ", mat, N);
#ifdef __amd64__
    test(LU_sse_fma, "SSE FMA: ", mat, N);
    test(LU_avx, "AVX: ", mat, N);
    test(LU_avx_aligned, "AVX Aligned: ", mat, N);
#endif
#ifdef __AVX512F__
    test(LU_avx512, "AVX512: ", mat, N);
    test(LU_avx512_aligned, "AVX512 Aligned: ", mat, N);
#endif
    cout << endl;
#else
    LU(mat, N);
    LU_simd(mat, N);
#ifdef __amd64__
    LU_avx(mat, N);
    LU_avx_aligned(mat, N);
#endif
#endif

    return 0;
}
#include <pthread.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <string.h>

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

#ifdef __amd64__
#include <immintrin.h>
#include "NEON_2_SSE.h"
#define USE_SSE4
#endif

#define ZERO (float)1e-5
#ifndef N
#define N 4096
#define NUM_THREADS 20
#endif
#define REPT 1
#define ele_t float
// #define DEBUG

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
    cout << time_used << ',';
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

struct LU_data
{
    int th;
    pthread_mutex_t finished = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t startNext = PTHREAD_MUTEX_INITIALIZER;
    ele_t (*mat)[N][N];
    int n;
    int i, begin, nLines; // 当前行、开始消去行、结束消去行
};

void *subthread_LU(void *_params)
{
    LU_data *params = (LU_data *)_params;
    int i = params->i;
    int n = params->n;
    float32x4_t mat_j, mat_i, div4;
    // cout << "mat_j addr: " << &mat_j << endl;
    // cout << "mat_i addr: " << &mat_i << endl;
    // cout << "div4 addr: " << &div4 << endl;
    for (int j = params->begin; j < params->begin + params->nLines; j++)
    {
        if ((*params->mat)[i][i] == 0)
            continue;
        ele_t div = (*params->mat)[j][i] / (*params->mat)[i][i];
        div4 = vmovq_n_f32(div);
        for (int k = i / 4 * 4; k < n; k += 4)
        {
            mat_j = vld1q_f32((*params->mat)[j] + k);
            mat_i = vld1q_f32((*params->mat)[i] + k);
            vst1q_f32((*params->mat)[j] + k, vmlsq_f32(mat_j, div4, mat_i));
        }
    }
}

// void *subthread_LU(void *_params)
// {
//     LU_data *params = (LU_data *)_params;
//     int i = params->i;
//     int n = params->n;
//     float(*mat)[N][N] = params->mat;
//     float i_mat[N];
//     // float n_mat[params->nLines][N];
//     // cout << params->begin << ' ' << params->nLines << endl;
//     for (int j = params->begin; j < params->begin + params->nLines; j++)
//     {
//         if ((*mat)[i][i] == 0)
//             continue;
//         ele_t div = (*mat)[j][i] / (*mat)[i][i];
//         // cout << new_mat[j][i] << '/' << new_mat[i][i] << '=' << div << endl;
//         for (int k = i; k < n; k++)
//             (*mat)[j][k] -= (*mat)[i][k] * div;
//     }
//     // cout << "i: " << params->i << endl;
//     // cout << "i: " << i << endl;
//     // cout << "begin: " << params->begin << endl;
//     // cout << params->th << " finished" << endl << endl;
// }

void LU_pthread(ele_t mat[N][N], int n)
{
    memcpy(new_mat, mat, sizeof(ele_t) * N * N);
    pthread_t threads[NUM_THREADS];
    LU_data attr[NUM_THREADS];

    for (int i = 0; i < n; i++)
    {
        // for (int j = i + 1; j < n; j++)
        // 从第i+1行开始遍历，步进为线程数
        int nLines = (n - i - 1) / NUM_THREADS;
        // cout << "-------------------" << endl;
        if (nLines > 31)
        {
            for (int th = 0; th < NUM_THREADS; th++)
            {
                attr[th].th = th;
                attr[th].mat = &new_mat;
                attr[th].n = n;
                attr[th].i = i;
                attr[th].nLines = nLines;
                attr[th].begin = i + 1 + th * nLines;
                // cout << "attr addr: " << &attr << endl;
                // cout << "attr->i: " << attr.i << endl;
                // cout << th << " creating" << endl;
                int err = pthread_create(&threads[th], NULL, subthread_LU, (void *)&attr[th]);

                // cout << "th" << th << endl;
                // int c;
                // for(int i=0; i<100000000; i++) c++;
                if (err)
                {
                    cout << "failed to create thread[" << th << "]" << endl;
                    exit(-1);
                }
                // pthread_join(threads[th], NULL);
            }

            // 算掉无法被整除的最后几行
            for (int j = i + 1 + NUM_THREADS * ((n - i - 1) / NUM_THREADS); j < n; j++)
            {
                if (new_mat[i][i] == 0)
                    continue;
                ele_t div = new_mat[j][i] / new_mat[i][i];
                for (int k = i; k < n; k++)
                    new_mat[j][k] -= new_mat[i][k] * div;
            }

            for (int th = 0; th < NUM_THREADS; th++)
                pthread_join(threads[th], NULL);
            // cout << "all finished" << endl << endl;
        }
        else
        {
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

void *subthread_static_LU(void *_params)
{
    LU_data *params = (LU_data *)_params;
    int i = params->i;
    int n = params->n;
    float32x4_t mat_j, mat_i, div4;
    // cout << "mat_j addr: " << &mat_j << endl;
    // cout << "mat_i addr: " << &mat_i << endl;
    // cout << "div4 addr: " << &div4 << endl;
    while (true)
    {
        pthread_mutex_lock(&(params->startNext));
        i = params->i;
        n = params->n;
        for (int j = params->begin; j < params->begin + params->nLines; j++)
        {
            if ((*params->mat)[i][i] == 0)
                continue;
            ele_t div = (*params->mat)[j][i] / (*params->mat)[i][i];
            div4 = vmovq_n_f32(div);
            for (int k = i / 4 * 4; k < n; k += 4)
            {
                mat_j = vld1q_f32((*params->mat)[j] + k);
                mat_i = vld1q_f32((*params->mat)[i] + k);
                vst1q_f32((*params->mat)[j] + k, vmlsq_f32(mat_j, div4, mat_i));
            }
        }
        pthread_mutex_unlock(&(params->finished));
    }
}

void LU_static_thread(ele_t mat[N][N], int n)
{
    memcpy(new_mat, mat, sizeof(ele_t) * N * N);
    pthread_t threads[NUM_THREADS];
    LU_data attr[NUM_THREADS];

    for (int th = 0; th < NUM_THREADS; th++)
    {
        pthread_mutex_lock(&(attr[th].startNext));
        pthread_mutex_lock(&(attr[th].finished));
        int err = pthread_create(&threads[th], NULL, subthread_static_LU, (void *)&attr[th]);
        if (err)
        {
            cout << "failed to create thread[" << th << "]" << endl;
            exit(-1);
        }
    }

    for (int i = 0; i < n; i++)
    {
        // for (int j = i + 1; j < n; j++)
        // 从第i+1行开始遍历，步进为线程数
        int nLines = (n - i - 1) / NUM_THREADS;
        // cout << "-------------------" << endl;

        for (int th = 0; th < NUM_THREADS; th++)
        {
            attr[th].th = th;
            attr[th].mat = &new_mat;
            attr[th].n = n;
            attr[th].i = i;
            attr[th].nLines = nLines;
            attr[th].begin = i + 1 + th * nLines;
            pthread_mutex_unlock(&(attr[th].startNext));
        }

        // 算掉无法被整除的最后几行
        for (int j = i + 1 + NUM_THREADS * ((n - i - 1) / NUM_THREADS); j < n; j++)
        {
            if (new_mat[i][i] == 0)
                continue;
            ele_t div = new_mat[j][i] / new_mat[i][i];
            for (int k = i; k < n; k++)
                new_mat[j][k] -= new_mat[i][k] * div;
        }

        for (int th = 0; th < NUM_THREADS; th++)
            pthread_mutex_lock(&(attr[th].finished));
        // cout << "all finished" << endl << endl;
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
    ifstream data("gauss.dat", ios::in | ios::binary);
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
    if (NUM_THREADS == 1)
        test(LU_simd, "NEON/SSE: ", mat, N);
    // test(LU_pthread, "pthread: ", mat, N);
    else
        test(LU_static_thread, "static thread: ", mat, N);
#else
    cout << endl;
    LU(mat, N);
    cout << endl
         << endl;
    // LU_simd(mat, N);
    // cout << endl
    //      << endl;
    LU_pthread(mat, N);
    cout << endl
         << endl;
    LU_static_thread(mat, N);
#endif
    return 0;
}

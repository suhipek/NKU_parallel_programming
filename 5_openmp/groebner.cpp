#include <pthread.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <string>
#include <string.h>
#include <omp.h>

#define PHILOSOPHY 能跑就行

#ifdef __ARM_NEON
#include <arm_neon.h>
#define simd_reg_t uint32x4_t
#endif

#ifdef __amd64__
#include <immintrin.h>
// #include "../NEON_2_SSE.h"
#define USE_SSE4
#ifdef __AVX512F__
#define simd_reg_t __m512i
#else
#define simd_reg_t __m256i
#endif

#endif

#define DEBUG

#ifndef NUM_THREADS
#define NUM_THREADS 18
#endif

// #ifndef DATA
// #define DATA "../Groebner/1_130_22_8/"
// #define COL 130
// #define ELE 22
// #define ROW 8
// #endif

#ifndef DATA
#define DATA "../Groebner/2_254_106_53/"
#define COL 254
#define ELE 106
#define ROW 53
#endif

// #ifndef DATA
// #define DATA "./Groebner/3_562_170_53/"
// #define COL 562
// #define ELE 170
// #define ROW 53
// #endif

// #ifndef DATA
// #define DATA "../Groebner/4_1011_539_263/"
// #define COL 1011
// #define ELE 539
// #define ROW 263
// #endif

// #ifndef DATA
// #define DATA "../Groebner/6_3799_2759_1953/"
// #define COL 3799
// #define ELE 2759
// #define ROW 1953
// #endif

// #ifndef DATA
// #define DATA "../Groebner/7_8399_6375_4535/"
// #define COL 8399
// #define ELE 6375
// #define ROW 4535
// #endif

// #ifndef DATA
// #define DATA "../Groebner/8_23075_18748_14325/"
// #define COL 23075
// #define ELE 18748
// #define ROW 14325
// #endif

// #ifndef DATA
// #define DATA "../Groebner/9_37960_29304_14921/"
// #define COL 37960
// #define ELE 29304
// #define ROW 14921
// #endif

#define mat_t unsigned int
#define mat_L 32
#define REPT 1

#ifndef SEPR
#define SEPR endl
#endif

#ifndef OPT_CLAUSE
#define NUM_BLOCKS 4
#define OPT_CLAUSE schedule(guided)
#endif

using namespace std;

mat_t ele[COL][COL / mat_L + 1] __attribute__((aligned(64))) = {0};
mat_t row[ROW][COL / mat_L + 1] __attribute__((aligned(64))) = {0};

// mat_t ele_tmp[COL][COL / mat_L + 1] __attribute__((aligned(64))) = {0};
// mat_t row_tmp[ROW][COL / mat_L + 1] __attribute__((aligned(64))) = {0};

void test(void (*func)(mat_t[COL][COL / mat_L + 1], mat_t[ROW][COL / mat_L + 1]), const char *msg)
{
    timespec start, end;
    double time_used = 0;
    // cout << "result: " << func(arr, len) << "    ";
    clock_gettime(CLOCK_REALTIME, &start);
    // for (int i = 0; i < REPT*(int)pow(2,(20-(int)(log2(len)))); i++)
    for (int i = 0; i < REPT; i++)
        func(ele, row);
    clock_gettime(CLOCK_REALTIME, &end);
    time_used += (end.tv_sec - start.tv_sec) * 1000;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000;
    cout << time_used << SEPR;
}

void groebner(mat_t ele_tmp[COL][COL / mat_L + 1], mat_t row_tmp[ROW][COL / mat_L + 1])
{
    // ele=消元子，row=被消元行
    // memcpy(ele_tmp, ele, sizeof(mat_t) * COL * (COL / mat_L + 1));
    // memcpy(row_tmp, row, sizeof(mat_t) * ROW * (COL / mat_L + 1));
    for (int i = 0; i < ROW; i++)
    { // 遍历被消元行
        for (int j = COL - 1; j >= 0; j--)
        { // 遍历列
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
            { // 当前位置有元素，需要消元
                if (ele_tmp[j][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                { //找到对应消元子
                    for (int p = COL / mat_L; p >= 0; p--)
                        row_tmp[i][p] ^= ele_tmp[j][p];
                }
                else
                { // 找不到对应消元子，升格当前消元行
                    memcpy(ele_tmp[j], row_tmp[i], (COL / mat_L + 1) * sizeof(mat_t));
                    break;
                }
            }
        }
    }

#ifdef DEBUG
    for (int i = 0; i < ROW; i++)
    {
        cout << i << ": ";
        for (int j = COL; j >= 0; j--)
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                cout << j << ' ';
        cout << endl;
    }
#endif
}

void groebner_new(mat_t ele_tmp[COL][COL / mat_L + 1], mat_t row_tmp[ROW][COL / mat_L + 1])
{
    // ele=消元子，row=被消元行
    // memcpy(ele_tmp, ele, sizeof(mat_t) * COL * (COL / mat_L + 1));
    // memcpy(row_tmp, row, sizeof(mat_t) * ROW * (COL / mat_L + 1));

    bool upgraded[ROW] = {0};

    for (int j = COL - 1; j >= 0; j--)
    { // 遍历消元子
        if (!(ele_tmp[j][j / mat_L] & ((mat_t)1 << (j % mat_L))))
        { // 如果存在对应消元子则进行消元
            for (int i = 0; i < ROW; i++)
            { // 遍历被消元行
                if (upgraded[i])
                    continue;
                if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                {
                    memcpy(ele_tmp[j], row_tmp[i], (COL / mat_L + 1) * sizeof(mat_t));
                    upgraded[i] = true;
                    break;
                }
            }
        }
        for (int i = 0; i < ROW; i++)
        { // 遍历被消元行
            if (upgraded[i])
                continue;
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
            { // 如果当前行需要消元
                for (int p = COL / mat_L; p >= 0; p--)
                    row_tmp[i][p] ^= ele_tmp[j][p];
            }
        }
    }

#ifdef DEBUG
    for (int i = 0; i < ROW; i++)
    {
        // cout << i << ": ";
        for (int j = COL; j >= 0; j--)
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                cout << j << ' ';
        cout << endl;
    }
#endif
}

void groebner_omp(mat_t ele_tmp[COL][COL / mat_L + 1], mat_t row_tmp[ROW][COL / mat_L + 1])
{
    // ele=消元子，row=被消元行
    // memcpy(ele_tmp, ele, sizeof(mat_t) * COL * (COL / mat_L + 1));
    // memcpy(row_tmp, row, sizeof(mat_t) * ROW * (COL / mat_L + 1));

    bool upgraded[ROW] = {0};

#pragma omp parallel num_threads(NUM_THREADS)
    for (int j = COL - 1; j >= 0; j--)
    { // 遍历消元子
#pragma omp master
        if (!(ele_tmp[j][j / mat_L] & ((mat_t)1 << (j % mat_L))))
        { // 如果不存在对应消元子则进行升格
            for (int i = 0; i < ROW; i++)
            { // 遍历被消元行
                if (upgraded[i])
                    continue;
                if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                {
                    // printf("线程%d: 不存在消元子%d, 将被消元行%d升格\n", omp_get_thread_num(), j, i);
                    memcpy(ele_tmp[j], row_tmp[i], (COL / mat_L + 1) * sizeof(mat_t));
                    upgraded[i] = true;
                    break;
                }
            }
        }
#pragma omp barrier
#pragma omp for OPT_CLAUSE
        for (int i = 0; i < ROW; i++)
        { // 遍历被消元行
            if (upgraded[i])
                continue;
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
            { // 如果当前行需要消元
#pragma omp simd
                for (int p = 0; p <= COL / mat_L; p++)
                    row_tmp[i][p] ^= ele_tmp[j][p];
            }
        }
    }

#ifdef DEBUG
    for (int i = 0; i < ROW; i++)
    {
        // cout << i << ": ";
        for (int j = COL; j >= 0; j--)
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                cout << j << ' ';
        cout << endl;
    }
#endif
}

int main()
{
    // cout << (string)DATA + (string) "1.txt" << endl;
    ifstream data_ele((string)DATA + (string) "1.txt", ios::in);
    int temp, header;
    string line;

    for (int i = 0; i < ELE; i++)
    {
        getline(data_ele, line);
        istringstream line_iss(line);
        line_iss >> header;
        ele[header][header / mat_L] += (mat_t)1 << (header % mat_L);
        while (line_iss >> temp)
            ele[header][temp / mat_L] += (mat_t)1 << (temp % mat_L);
    }
    data_ele.close();

    ifstream data_row((string)DATA + (string) "2.txt", ios::in);
    for (int i = 0; i < ROW; i++)
    {
        getline(data_row, line);
        istringstream line_iss(line);
        while (line_iss >> temp)
            row[i][temp / mat_L] += (mat_t)1 << (temp % mat_L);
    }
    data_row.close();

#ifdef DEBUG
    // groebner(ele, row);
    // cout << endl
    //      << "end" << endl;
    // groebner_new(ele, row);
    // cout << endl
    //      << "end" << endl;
    groebner_omp(ele, row);
#else
    if (NUM_THREADS == 1)
        test(groebner, "common");
    else
        test(groebner_omp, "common");
#endif
    return 0;
}
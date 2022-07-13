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
#include "../xsimd/include/xsimd/xsimd.hpp"

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

// #define DEBUG

#ifndef NUM_THREADS
#define NUM_THREADS 1
#endif

// #ifndef DATA
// #define DATA "../Groebner/1_130_22_8/"
// #define COL 130
// #define ELE 22
// #define ROW 8
// #endif

// #ifndef DATA
// #define DATA "../Groebner/2_254_106_53/"
// #define COL 254
// #define ELE 106
// #define ROW 53
// #endif

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
// #define DATA "../Groebner/5_2362_1226_453/"
// #define COL 2362
// #define ELE 1226
// #define ROW 453
// #endif

#ifndef DATA
#define DATA "../Groebner/6_3799_2759_1953/"
#define COL 3799
#define ELE 2759
#define ROW 1953
#endif

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

// #ifndef DATA
// #define DATA "../Groebner/11_85401_5724_756/"
// #define COL 85401
// #define ELE 5724
// #define ROW 756
// #endif

#define mat_t unsigned int
#define mat_L 32
#define REPT 1

#ifndef SEPR
#define SEPR endl
#endif

#ifndef OPT_CLAUSE
#define NUM_BLOCKS 4
#define OPT_CLAUSE schedule(static)
#endif

using namespace std;

mat_t ele[COL][(COL / mat_L + 1) / 16 * 16 + 16] __attribute__((aligned(128))) = {0};
mat_t row[ROW][(COL / mat_L + 1) / 16 * 16 + 16] __attribute__((aligned(128))) = {0};

// mat_t ele_tmp[COL][(COL / mat_L + 1) / 16 * 16 + 16] __attribute__((aligned(64))) = {0};
// mat_t row_tmp[ROW][(COL / mat_L + 1) / 16 * 16 + 16] __attribute__((aligned(64))) = {0};

void test(void (*func)(mat_t[COL][(COL / mat_L + 1) / 16 * 16 + 16], mat_t[ROW][(COL / mat_L + 1) / 16 * 16 + 16]), const char *msg)
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

void groebner_new(mat_t ele_tmp[COL][(COL / mat_L + 1) / 16 * 16 + 16], mat_t row_tmp[ROW][(COL / mat_L + 1) / 16 * 16 + 16])
{
    // ele=消元子，row=被消元行
    // memcpy(ele_tmp, ele, sizeof(mat_t) * COL * ((COL / mat_L + 1) / 16 * 16 + 16));
    // memcpy(row_tmp, row, sizeof(mat_t) * ROW * ((COL / mat_L + 1) / 16 * 16 + 16));

    bool upgraded[ROW] = {0};

    for (int j = COL - 1; j >= 0; j--)
    { // 遍历消元子
        if (!(ele_tmp[j][j / mat_L] & ((mat_t)1 << (j % mat_L))))
        { // 如果不存在对应消元子则进行升格
            for (int i = 0; i < ROW; i++)
            { // 遍历被消元行
                if (upgraded[i])
                    continue;
                if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                {
                    memcpy(ele_tmp[j], row_tmp[i], ((COL / mat_L + 1) / 16 * 16 + 16) * sizeof(mat_t));
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

template<class Arch, class Mode>
void groebner_new_simd(mat_t ele_tmp[COL][(COL / mat_L + 1) / 16 * 16 + 16], mat_t row_tmp[ROW][(COL / mat_L + 1) / 16 * 16 + 16])
{
    // ele=消元子，row=被消元行
    // memcpy(ele_tmp, ele, sizeof(mat_t) * COL * ((COL / mat_L + 1) / 16 * 16 + 16));
    // memcpy(row_tmp, row, sizeof(mat_t) * ROW * ((COL / mat_L + 1) / 16 * 16 + 16));

    bool upgraded[ROW] = {0};
    using mat_simd_t = xsimd::batch<unsigned int, Arch>;
    std::size_t inc = mat_simd_t::size;
    mat_simd_t ele_vec, row_vec;

    for (int j = COL - 1; j >= 0; j--)
    { // 遍历消元子
        if (!(ele_tmp[j][j / mat_L] & ((mat_t)1 << (j % mat_L))))
        { // 如果不存在对应消元子则进行升格
            for (int i = 0; i < ROW; i++)
            { // 遍历被消元行
                if (upgraded[i])
                    continue;
                if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                {
                    memcpy(ele_tmp[j], row_tmp[i], ((COL / mat_L + 1) / 16 * 16 + 16) * sizeof(mat_t));
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
                int p=0;
                for(; p <= COL / mat_L; p+=inc)
                {
                    ele_vec = mat_simd_t::load_unaligned(&ele_tmp[j][p]);
                    row_vec = mat_simd_t::load_unaligned(&row_tmp[i][p]);
                    row_vec ^= ele_vec;
                    xsimd::store_unaligned(&row_tmp[i][p], row_vec);
                }
                for(; p <= COL / mat_L; p++)
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
    //groebner_new(ele, row);
    // cout << endl
    //      << "end" << endl;
    // groebner_omp(ele, row);
    groebner_new_simd<xsimd::avx, xsimd::aligned_mode>(ele, row);
#else
    // test(groebner_new, "common");
    // test(groebner_new_simd<xsimd::sse3, xsimd::unaligned_mode>, "avx");
    // test(groebner_new_simd<xsimd::avx, xsimd::unaligned_mode>, "avx");
    test(groebner_new_simd<xsimd::avx2, xsimd::aligned_mode>, "avx");
#endif
    return 0;
}
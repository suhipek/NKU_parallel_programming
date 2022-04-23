#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <string>
#include <string.h>

#define PHILOSOPHY 能跑就行

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

#ifdef __amd64__
#include <immintrin.h>
#include "NEON_2_SSE.h"
#define USE_SSE4
#endif

// #ifndef DATA
// #define DATA "./Groebner/1_130_22_8/"
// #define COL 130
// #define ELE 22
// #define ROW 8
// #endif

// #ifndef DATA
// #define DATA "./Groebner/2_254_106_53/"
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
// #define DATA "./Groebner/4_1011_539_263/"
// #define COL 1011
// #define ELE 539
// #define ROW 263
// #endif

#ifndef DATA
#define DATA "./Groebner/6_3799_2759_1953/"
#define COL 3799
#define ELE 2759
#define ROW 1953
#endif

#define mat_t unsigned int
#define mat_L 32
#define REPT 4

using namespace std;

mat_t ele[COL][COL / mat_L + 1] = {0};
mat_t row[ROW][COL / mat_L + 1] = {0};

#ifndef ALIGN
mat_t ele_tmp[COL][COL / mat_L + 1] = {0};
mat_t row_tmp[ROW][COL / mat_L + 1] = {0};
#else
mat_t ele_tmp[COL][(COL / mat_L + 1) / 16 * 16 + 16] __attribute__((aligned(64))) = {0};
mat_t row_tmp[ROW][(COL / mat_L + 1) / 16 * 16 + 16] __attribute__((aligned(64))) = {0};
#endif

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
    time_used += end.tv_sec - start.tv_sec;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
    cout << time_used << ',';
}

void groebner(mat_t ele[COL][COL / mat_L + 1], mat_t row[ROW][COL / mat_L + 1])
{
    // ele=消元子，row=被消元行
    memcpy(ele_tmp, ele, sizeof(mat_t) * COL * (COL / mat_L + 1));
    memcpy(row_tmp, row, sizeof(mat_t) * ROW * (COL / mat_L + 1));
    for (int i = 0; i < ROW; i++)
    {
        for (int j = COL; j >= 0; j--)
        {
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
            {
                if (ele_tmp[j][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                {
                    for (int p = COL / mat_L; p >= 0; p--)
                        row_tmp[i][p] ^= ele_tmp[j][p];
                }
                else
                {
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

void groebner_simd(mat_t ele[COL][COL / mat_L + 1], mat_t row[ROW][COL / mat_L + 1])
{
    // ele=消元子，row=被消元行
    memcpy(ele_tmp, ele, sizeof(mat_t) * COL * (COL / mat_L + 1));
    memcpy(row_tmp, row, sizeof(mat_t) * ROW * (COL / mat_L + 1));
    uint32x4_t row_i, ele_j;
    for (int i = 0; i < ROW; i++)
    {
        for (int j = COL; j >= 0; j--)
        {
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
            {
                if (ele_tmp[j][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                {
                    for (int p = 0; p < COL / 128; p++)
                    {
                        row_i = vld1q_u32(row_tmp[i] + p * 4);
                        ele_j = vld1q_u32(ele_tmp[j] + p * 4);
                        vst1q_u32(row_tmp[i] + p * 4, veorq_u32(row_i, ele_j));
                    }
                    for (int k = COL / 128 * 4; k <= COL / mat_L; k++)
                        row_tmp[i][k] ^= ele_tmp[j][k];
                }
                else
                {
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

#ifdef __amd64__
void groebner_avx(mat_t ele[COL][COL / mat_L + 1], mat_t row[ROW][COL / mat_L + 1])
{
    // ele=消元子，row=被消元行
    memcpy(ele_tmp, ele, sizeof(mat_t) * COL * (COL / mat_L + 1));
    memcpy(row_tmp, row, sizeof(mat_t) * ROW * (COL / mat_L + 1));
    __m256i row_i, ele_j;
    for (int i = 0; i < ROW; i++)
    {
        for (int j = COL; j >= 0; j--)
        {
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
            {
                if (ele_tmp[j][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                {
                    for (int p = 0; p < COL / 256; p++)
                    {
                        #ifdef ALIGN
                        row_i = _mm256_load_si256((__m256i *)(row_tmp[i] + p * 8));
                        ele_j = _mm256_load_si256((__m256i *)(ele_tmp[j] + p * 8));
                        _mm256_store_si256((__m256i *)(row_tmp[i] + p * 8), _mm256_xor_si256(row_i, ele_j));
                        #else
                        row_i = _mm256_loadu_si256((__m256i *)(row_tmp[i] + p * 8));
                        ele_j = _mm256_loadu_si256((__m256i *)(ele_tmp[j] + p * 8));
                        _mm256_storeu_si256((__m256i *)(row_tmp[i] + p * 8), _mm256_xor_si256(row_i, ele_j));
                        #endif
                    }
                    for (int k = COL / 256 * 8; k <= COL / mat_L; k++)
                        row_tmp[i][k] ^= ele_tmp[j][k];
                }
                else
                {
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
#endif

#ifdef __AVX512F__
void groebner_avx512(mat_t ele[COL][COL / mat_L + 1], mat_t row[ROW][COL / mat_L + 1])
{
    // ele=消元子，row=被消元行
    memcpy(ele_tmp, ele, sizeof(mat_t) * COL * (COL / mat_L + 1));
    memcpy(row_tmp, row, sizeof(mat_t) * ROW * (COL / mat_L + 1));
    __m512i row_i, ele_j;
    for (int i = 0; i < ROW; i++)
    {
        for (int j = COL; j >= 0; j--)
        {
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
            {
                if (ele_tmp[j][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                {
                    for (int p = 0; p < COL / 512; p++)
                    {
                        #ifdef ALIGN
                        row_i = _mm512_load_si512((__m512i *)(row_tmp[i] + p * 16));
                        ele_j = _mm512_load_si512((__m512i *)(ele_tmp[j] + p * 16));
                        _mm512_store_si512((__m512i *)(row_tmp[i] + p * 16), _mm512_xor_si512(row_i, ele_j));
                        #else
                        row_i = _mm512_loadu_si512((__m512i *)(row_tmp[i] + p * 16));
                        ele_j = _mm512_loadu_si512((__m512i *)(ele_tmp[j] + p * 16));
                        _mm512_storeu_si512((__m512i *)(row_tmp[i] + p * 16), _mm512_xor_si512(row_i, ele_j));
                        #endif
                    }
                    for (int k = COL / 512 * 16; k <= COL / mat_L; k++)
                        row_tmp[i][k] ^= ele_tmp[j][k];
                }
                else
                {
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
#endif

int main()
{
    //cout << (string)DATA + (string) "1.txt" << endl;
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
    groebner(ele, row);
    cout << endl << "end" << endl;
    groebner_simd(ele, row);
    cout << endl << "end" << endl;
    groebner_avx(ele, row);
#else
    test(groebner, "common");
    test(groebner_simd, "simd");
#ifdef __amd64__
    test(groebner_avx, "avx");
#endif
#ifdef __AVX512F__
    test(groebner_avx512, "avx512");
#endif
#endif

    // groebner(ele, row);
    //  for (int i = 0; i < COL; i++)
    //  {
    //      cout << i << ": ";
    //      for (int j = COL; j >= 0; j--)
    //          if (ele[i][j / mat_L] & (mat_t)(1 << (j % mat_L)))
    //              cout << j << ' ';
    //      cout << endl;
    //  }
    return 0;
}
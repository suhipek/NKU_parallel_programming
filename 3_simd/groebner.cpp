#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <string>
#include <string.h>

#ifdef __arm__
#include <arm_neon.h>
#endif

// #ifdef __amd64__
// #include <immintrin.h>
// #include "NEON_2_SSE.h"
// #define USE_SSE4
// #endif

#define DATA "./Groebner/1_130_22_8/"
#define COL 130
#define ELE 22
#define ROW 8
#define mat_t unsigned int
#define mat_L 32
using namespace std;

mat_t ele[COL][COL / mat_L + 1] = {0};
mat_t row[ROW][COL / mat_L + 1] = {0};

mat_t ele_tmp[COL][COL / mat_L + 1] = {0};
mat_t row_tmp[ROW][COL / mat_L + 1] = {0};

// void test(void (*func)(ele_t[N][N], int), const char *msg, ele_t mat[N][N], int len)
// {
//     timespec start, end;
//     double time_used = 0;
//     // cout << "result: " << func(arr, len) << "    ";
//     clock_gettime(CLOCK_REALTIME, &start);
//     // for (int i = 0; i < REPT*(int)pow(2,(20-(int)(log2(len)))); i++)
//     for (int i = 0; i < REPT; i++)
//         func(mat, len);
//     clock_gettime(CLOCK_REALTIME, &end);
//     time_used += end.tv_sec - start.tv_sec;
//     time_used += double(end.tv_nsec - start.tv_nsec) / 1000000000;
//     cout << time_used << ',';
// }

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
                        row_tmp[i][p] = (row_tmp[i][p] ^ ele_tmp[j][p]);
                }
                else
                {
                    memcpy(ele_tmp[j], row_tmp[i], (COL / mat_L + 1) * sizeof(mat_t));
                    break;
                }
            }
        }
    }

    for (int i = 0; i < ROW; i++)
    {
        for (int j = COL; j >= 0; j--)
            if (row_tmp[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                cout << j << ' ';
        cout << endl;
    }
}

int main()
{
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

    groebner(ele, row);
    // for (int i = 0; i < COL; i++)
    // {
    //     cout << i << ": ";
    //     for (int j = COL; j >= 0; j--)
    //         if (ele[i][j / mat_L] & (mat_t)(1 << (j % mat_L)))
    //             cout << j << ' ';
    //     cout << endl;
    // }
    return 0;
}
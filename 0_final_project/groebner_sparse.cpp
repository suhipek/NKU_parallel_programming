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
#include <vector>
#include <unordered_map>

using namespace std;

#define mat_t int

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

#ifndef DATA
#define DATA "../Groebner/5_2362_1226_453/"
#define COL 2362
#define ELE 1226
#define ROW 453
#endif

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
// #define DATA "../Groebner/11_85401_5724_756/"
// #define COL 85401
// #define ELE 5724
// #define ROW 756
// #endif

// #define DEBUG

void groebner_sparse(array<vector<mat_t>, COL> ele, array<vector<mat_t>, ROW> row)
{
    bool upgraded[ROW] = {0};
    mat_t buffer[COL];

    for (mat_t j = COL - 1; j >= 0; j--)
    { // 遍历消元子
        if (ele[j].size() == 0)
        { // 如果不存在对应消元子则进行升格
            for (mat_t i = 0; i < ROW; i++)
            { // 遍历被消元行
                if (upgraded[i])
                    continue;
                if (row[i][0] == j)
                {
                    ele[j] = row[i];
                    upgraded[i] = true;
                    break;
                }
            }
        }
        for (mat_t i = 0; i < ROW; i++)
        { // 遍历被消元行
            if (upgraded[i])
                continue;
            if (row[i][0] == j)
            { // 如果当前行需要消元
                mat_t pRow = 0, pEle = 0;
                int rowMax = row[i].size();
                int eleMax = ele[j].size();
                mat_t row_i_p, ele_j_p;
                vector<mat_t> result;

                mat_t last = 0;

                while (pRow < rowMax && pEle < eleMax)
                {
                    row_i_p = row[i][pRow];
                    ele_j_p = ele[j][pEle];
                    
                    pRow += row_i_p >= ele_j_p ? 1 : 0;
                    pEle += ele_j_p >= row_i_p ? 1 : 0;
                    buffer[last] = row_i_p > ele_j_p ? row_i_p:ele_j_p;
                    // buffer[last] = row_i_p;
                    //last++;
                    last += row_i_p == ele_j_p ? 0 : 1;
                }

                for (; pRow < rowMax; pRow++)
                {
                    buffer[last]=row[i][pRow];
                    last++;
                }
                for (; pEle < eleMax; pEle++)
                {
                    buffer[last]=ele[j][pEle];
                    last++;
                }
                row[i].resize(last);
                memcpy(&(row[i][0]), buffer, last*sizeof(mat_t));
            }
        }
    }
#ifdef DEBUG
    for (auto i : row)
    {
        // cout << i[0] << ": ";
        for (auto j : i)
            cout << j << ' ';
        cout << endl;
    }
#endif
}

int main()
{
    ifstream data_ele((string)DATA + (string) "1.txt", ios::in);
    mat_t temp, header;
    string line;
    array<vector<mat_t>, COL> ele;
    array<vector<mat_t>, ROW> row;

    for (mat_t i = 0; i < ELE; i++)
    {
        getline(data_ele, line);
        istringstream line_iss(line);
        line_iss >> header;
        ele[header].push_back(header);
        while (line_iss >> temp)
            ele[header].push_back(temp);
    }
    data_ele.close();

    ifstream data_row((string)DATA + (string) "2.txt", ios::in);
    for (mat_t i = 0; i < ROW; i++)
    {
        getline(data_row, line);
        istringstream line_iss(line);
        while (line_iss >> temp)
            row[i].push_back(temp);
    }
    data_row.close();

    timespec start, end;
    double time_used = 0;
    clock_gettime(CLOCK_REALTIME, &start);

    groebner_sparse(ele, row);

    clock_gettime(CLOCK_REALTIME, &end);
    time_used += (end.tv_sec - start.tv_sec) * 1000;
    time_used += double(end.tv_nsec - start.tv_nsec) / 1000000;
    cout << time_used << endl;

    // for(auto i : ele)
    // {
    //     // cout << i.first << ": ";
    //     for(auto j : i.second)
    //         cout << j << ' ';
    //     cout << endl;
    // }

    // cout << "end" << endl;

    // for(auto i : row)
    // {
    //     // cout << i[0] << ": ";
    //     for(auto j : i)
    //         cout << j << ' ';
    //     cout << endl;
    // }

    // int row[9] = {123, 100, 99, 32, 18, 10, 5, 2, -1};
    // int ele[8] = {123, 101, 90, 24, 18, 6, 1, -1};
    // vector<int> result;
    // int pRow = 0, pEle = 0;

    // while(pRow < 9 && pEle < 8)
    // {
    //     if(row[pRow] == ele[pEle])
    //     {
    //         pRow++;
    //         pEle++;
    //     }
    //     else if(row[pRow] < ele[pEle])
    //     {
    //         result.push_back(ele[pEle]);
    //         pEle++;
    //     }
    //     else
    //     {
    //         result.push_back(row[pRow]);
    //         pRow++;
    //     }
    // }
    // for(int i = 0; i < result.size(); i++)
    // {
    //     cout << result[i] << " ";
    // }
    return 0;
}
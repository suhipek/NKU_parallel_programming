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
#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

namespace mpi = boost::mpi;
using namespace std;

#define mat_t int

// #ifndef DATA
// #define DATA "/home/suhipek/NKU_parallel_programming/Groebner/1_130_22_8/"
// #define COL 130
// #define ELE 22
// #define ROW 8
// #endif

// #ifndef DATA
// #define DATA "/home/suhipek/NKU_parallel_programming/Groebner/2_254_106_53/"
// #define COL 254
// #define ELE 106
// #define ROW 53
// #endif

// #ifndef DATA
// #define DATA "/home/suhipek/NKU_parallel_programming/Groebner/5_2362_1226_453/"
// #define COL 2362
// #define ELE 1226
// #define ROW 453
// #endif

#ifndef DATA
#define DATA "/home/suhipek/NKU_parallel_programming/Groebner/6_3799_2759_1953/"
#define COL 3799
#define ELE 2759
#define ROW 1953
#endif

// #ifndef DATA
// #define DATA "/home/suhipek/NKU_parallel_programming/Groebner/7_8399_6375_4535/"
// #define COL 8399
// #define ELE 6375
// #define ROW 4535
// #endif

// #ifndef DATA
// #define DATA "/home/suhipek/NKU_parallel_programming/Groebner/11_85401_5724_756/"
// #define COL 85401
// #define ELE 5724
// #define ROW 756
// #endif

// #define DEBUG

mpi::communicator world;

void run_master(array<vector<mat_t>, COL> ele, array<vector<mat_t>, ROW> row)
{
    bool upgraded[ROW] = {0};
    mat_t buffer[COL];
    int n_workload = ROW / world.size() + 1;
    for (auto i_ele : ele)
    {
        broadcast(world, i_ele, 0);
    }
    for (int i = 0; i / n_workload + 1 < world.size(); i++)
    {
        int th = i / n_workload + 1;
        world.send(th, 0, row[i]);
    }
    for (mat_t j = COL - 1; j >= 0; j--)
    { // 遍历消元子
        if (ele[j].size() == 0)
        { // 如果不存在对应消元子则进行升格
            int tobe_upgraded = (1 << 31) - 1, real_upgraded = 0;
            for (mat_t i = n_workload * (world.size() - 1); i < ROW; i++)
            { // 遍历被消元行
                if (upgraded[i])
                    continue;
                if (row[i][0] == j)
                {
                    // ele[j] = row[i];
                    // upgraded[i] = true;
                    tobe_upgraded = i;
                    break;
                }
            }
            mpi::reduce(world, tobe_upgraded, real_upgraded, mpi::minimum<mat_t>(), 0);
            mpi::broadcast(world, real_upgraded, 0);
            if (real_upgraded == (1 << 31) - 1)
                continue;
            int th_with_min = (real_upgraded / n_workload + 1) % world.size();
            if (tobe_upgraded != (1 << 31) - 1)
                ele[j] = row[tobe_upgraded];
            mpi::broadcast(world, ele[j], th_with_min);
        }
        world.barrier();
        for (mat_t i = n_workload * (world.size() - 1); i < ROW; i++)
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
                    buffer[last] = row_i_p > ele_j_p ? row_i_p : ele_j_p;
                    last += row_i_p == ele_j_p ? 0 : 1;
                }

                for (; pRow < rowMax; pRow++)
                {
                    buffer[last] = row[i][pRow];
                    last++;
                }
                for (; pEle < eleMax; pEle++)
                {
                    buffer[last] = ele[j][pEle];
                    last++;
                }
                row[i].resize(last);
                std::memcpy(&(row[i][0]), buffer, last * sizeof(mat_t));
            }
        }
        world.barrier();
    }
    for (int i = 0; i / n_workload + 1 < world.size(); i++)
    {
        int th = i / n_workload + 1;
        world.recv(th, 0, row[i]);
    }

#ifdef DEBUG
    for (auto i : row)
    {
        // std::cout << i[0] << ": ";
        for (auto j : i)
            std::cout << j << ' ';
        std::cout << endl;
    }
#endif
}
void run_slave()
{
    bool upgraded[ROW] = {0};
    mat_t buffer[COL];
    int n_workload = ROW / world.size() + 1;

    array<vector<mat_t>, COL> ele;
    vector<vector<mat_t>> row;
    row.resize(n_workload);
    for (int i=0; i<COL; i++)
    {
        broadcast(world, ele[i], 0);
    }
    int iBegin = (world.rank() - 1) * n_workload;
    int iEnd = iBegin + n_workload;
    for (int i = iBegin; i < iEnd; i++)
    {
        world.recv(0, 0, row[i - iBegin]);
    }

    for (mat_t j = COL - 1; j >= 0; j--)
    { // 遍历消元子
        if (ele[j].size() == 0)
        { // 如果不存在对应消元子则进行升格
            int tobe_upgraded = (1 << 31) - 1, real_upgraded = 0;
            for (int i = 0; i < n_workload; i++)
            { // 遍历被消元行
                if (upgraded[i])
                    continue;
                if (row[i][0] == j)
                {
                    // ele[j] = row[i];
                    // upgraded[i] = true;
                    tobe_upgraded = i + iBegin;
                    break;
                }
            }
            mpi::reduce(world, tobe_upgraded, real_upgraded, mpi::minimum<mat_t>(), 0);
            mpi::broadcast(world, real_upgraded, 0);
            if (real_upgraded == (1 << 31) - 1)
                continue;
            int th_with_min = (real_upgraded / n_workload + 1) % world.size();
            if (tobe_upgraded != (1 << 31) - 1)
                ele[j] = row[tobe_upgraded - iBegin];
            mpi::broadcast(world, ele[j], th_with_min);
        }
        world.barrier();
        for (mat_t i = iBegin; i < iEnd; i++)
        { // 遍历被消元行
            if (upgraded[i])
                continue;
            if (row[i - iBegin][0] == j)
            { // 如果当前行需要消元
                mat_t pRow = 0, pEle = 0;
                int rowMax = row[i - iBegin].size();
                int eleMax = ele[j].size();
                mat_t row_i_p, ele_j_p;
                vector<mat_t> result;

                mat_t last = 0;

                while (pRow < rowMax && pEle < eleMax)
                {
                    row_i_p = row[i - iBegin][pRow];
                    ele_j_p = ele[j][pEle];
                    pRow += row_i_p >= ele_j_p ? 1 : 0;
                    pEle += ele_j_p >= row_i_p ? 1 : 0;
                    buffer[last] = row_i_p > ele_j_p ? row_i_p : ele_j_p;
                    last += row_i_p == ele_j_p ? 0 : 1;
                }

                for (; pRow < rowMax; pRow++)
                {
                    buffer[last] = row[i - iBegin][pRow];
                    last++;
                }
                for (; pEle < eleMax; pEle++)
                {
                    buffer[last] = ele[j][pEle];
                    last++;
                }
                row[i - iBegin].resize(last);
                std::memcpy(&(row[i - iBegin][0]), buffer, last * sizeof(mat_t));
            }
        }
        world.barrier();
    }
    for (int i = iBegin; i < iEnd; i++)
    {
        world.send(0, 0, row[i - iBegin]);
    }
}

int main()
{
    mpi::environment env;

    if (world.rank() == 0)
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

        run_master(ele, row);

        clock_gettime(CLOCK_REALTIME, &end);
        time_used += (end.tv_sec - start.tv_sec) * 1000;
        time_used += double(end.tv_nsec - start.tv_nsec) / 1000000;
        std::cout << time_used << endl;
    }
    else
    {
        run_slave();
    }
    return 0;
}
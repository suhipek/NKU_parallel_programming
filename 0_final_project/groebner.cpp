#include <mpi.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <string.h>
#include <omp.h>
#include <unistd.h>
#include <boost/mpi.hpp>
#include "../xsimd/include/xsimd/xsimd.hpp"

using namespace std;

// #define DEBUG

#ifndef NUM_THREADS
#define NUM_THREADS 16
#endif

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
// #define DATA "/home/suhipek/NKU_parallel_programming/Groebner/3_562_170_53/"
// #define COL 562
// #define ELE 170
// #define ROW 53
// #endif

// #ifndef DATA
// #define DATA "/home/suhipek/NKU_parallel_programming/Groebner/4_1011_539_263/"
// #define COL 1011
// #define ELE 539
// #define ROW 263
// #endif

// #ifndef DATA
// #define DATA "/home/suhipek/NKU_parallel_programming/Groebner/6_3799_2759_1953/"
// #define COL 3799
// #define ELE 2759
// #define ROW 1953
// #endif

#ifndef DATA
#define DATA "/home/suhipek/NKU_parallel_programming/Groebner/7_8399_6375_4535/"
#define COL 8399
#define ELE 6375
#define ROW 4535
#endif

// #ifndef DATA
// #define DATA "/home/suhipek/NKU_parallel_programming/Groebner/8_23075_18748_14325/"
// #define COL 23075
// #define ELE 18748
// #define ROW 14325
// #endif

#define mat_t unsigned int
#define mat_L 32
#define REPT 1
#define LEN_LINE ((COL / mat_L + 1) / 16 * 16 + 16)
#define SIMD_INST_SET xsimd::avx2
#define SIMD_ALIGN xsimd::aligned_mode()
// #define DEBUG

int world_size, world_rank;

void run_master(mat_t (*ele)[LEN_LINE], mat_t (*row)[LEN_LINE])
{
    bool upgraded[ROW] = {0};

    using mat_simd_t = xsimd::batch<unsigned int, SIMD_INST_SET>;
    std::size_t simd_inc = mat_simd_t::size;
    mat_simd_t ele_vec, row_vec;

    int n_workload = ROW / world_size + 1;
    for (int th = 1; th < world_size; th++)
    {
        int i_offset = n_workload * (th - 1);
        MPI_Send(row + i_offset, n_workload * (LEN_LINE) * sizeof(mat_t),
                 MPI_BYTE, th, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(ele, COL * (LEN_LINE) * sizeof(mat_t),
              MPI_BYTE, 0, MPI_COMM_WORLD);

    for (int j = COL - 1; j >= 0; j--)
    { // 遍历消元子
        if (!(ele[j][j / mat_L] & ((mat_t)1 << (j % mat_L))))
        { // 如果不存在对应消元子则进行升格
            int tobe_upgraded = (1 << 31) - 1, real_upgraded = 0;
            for (int i = n_workload * (world_size - 1); i < ROW; i++)
            { // 遍历被消元行
                if (upgraded[i])
                    continue;
                if (row[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                {
                    tobe_upgraded = i;
                    break;
                }
            }
            MPI_Reduce(&tobe_upgraded, &real_upgraded, 1, MPI_INT,
                       MPI_MIN, 0, MPI_COMM_WORLD);
            MPI_Bcast(&real_upgraded, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (real_upgraded == (1 << 31) - 1)
                continue;
            int th_with_min = (real_upgraded / n_workload + 1) % world_size;
            if (tobe_upgraded != (1 << 31) - 1)
                memcpy(ele[j], row[tobe_upgraded],
                       (LEN_LINE) * sizeof(mat_t));
            MPI_Bcast(ele[j], (LEN_LINE) * sizeof(mat_t),
                      MPI_BYTE, th_with_min, MPI_COMM_WORLD);
            upgraded[real_upgraded] = true;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // #pragma omp parallel for num_threads(4)
        for (int i = n_workload * (world_size - 1); i < ROW; i++)
        { // 遍历被消元行
            if (upgraded[i])
                continue;
            if (row[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
            { // 如果当前行需要消元
                int p = 0;
                for (; p <= COL / mat_L; p += simd_inc)
                {
                    ele_vec = mat_simd_t::load(&ele[j][p], SIMD_ALIGN);
                    row_vec = mat_simd_t::load(&row[i][p], SIMD_ALIGN);
                    row_vec ^= ele_vec;
                    xsimd::store(&row[i][p], row_vec, SIMD_ALIGN);
                }
                for (; p <= COL / mat_L; p++)
                    row[i][p] ^= ele[j][p];
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    for (int th = 1; th < world_size; th++)
    { // 回收被消元行
        int i_offset = n_workload * (th - 1);
        MPI_Recv(row + i_offset, n_workload * (LEN_LINE) * sizeof(mat_t),
                 MPI_BYTE, th, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

#ifdef DEBUG
    for (int i = 0; i < ROW; i++)
    {
        for (int j = COL; j >= 0; j--)
            if (row[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                cout << j << ' ';
        cout << endl;
    }
#endif
}

void run_slave()
{
    bool upgraded[ROW] = {0};
    int n_workload = ROW / world_size + 1;
    int i_offset = n_workload * (world_rank - 1);

    using mat_simd_t = xsimd::batch<unsigned int, SIMD_INST_SET>;
    std::size_t simd_inc = mat_simd_t::size;
    mat_simd_t ele_vec, row_vec;

    mat_t(*ele)[LEN_LINE] =
        (mat_t(*)[LEN_LINE])aligned_alloc(128, sizeof(mat_t) * COL * (LEN_LINE));
    mat_t(*row)[LEN_LINE] =
        (mat_t(*)[LEN_LINE])aligned_alloc(128, n_workload * (LEN_LINE) * sizeof(mat_t));
    MPI_Recv(row, n_workload * (LEN_LINE) * sizeof(mat_t),
             MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Bcast(ele, COL * (LEN_LINE) * sizeof(mat_t),
              MPI_BYTE, 0, MPI_COMM_WORLD);

    for (int j = COL - 1; j >= 0; j--)
    { // 遍历消元子
        if (!(ele[j][j / mat_L] & ((mat_t)1 << (j % mat_L))))
        { // 如果不存在对应消元子则进行升格
            int tobe_upgraded = (1 << 31) - 1, real_upgraded = 0;
            for (int i = 0; i < n_workload; i++)
            { // 遍历被消元行
                if (upgraded[i + i_offset])
                    continue;
                if (row[i][j / mat_L] & ((mat_t)1 << (j % mat_L)))
                {
                    tobe_upgraded = i + i_offset;
                    // memcpy(ele[j], row[i], (LEN_LINE) * sizeof(mat_t));
                    // upgraded[i] = true;
                    break;
                }
            }
            MPI_Reduce(&tobe_upgraded, &real_upgraded, 1, MPI_INT,
                       MPI_MIN, 0, MPI_COMM_WORLD);
            MPI_Bcast(&real_upgraded, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (real_upgraded == (1 << 31) - 1)
                continue;
            int th_with_min = (real_upgraded / n_workload + 1) % world_size;
            if (tobe_upgraded != (1 << 31) - 1)
                memcpy(ele[j], row[tobe_upgraded - i_offset],
                       (LEN_LINE) * sizeof(mat_t));
            MPI_Bcast(ele[j], (LEN_LINE) * sizeof(mat_t),
                      MPI_BYTE, th_with_min, MPI_COMM_WORLD);
            upgraded[real_upgraded] = true;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // #pragma omp parallel for num_threads(4)
        for (int i = i_offset; i < i_offset + n_workload; i++)
        {
            if (upgraded[i])
                continue;
            if (row[i - i_offset][j / mat_L] & ((mat_t)1 << (j % mat_L)))
            { // 如果当前行需要消元
                int p = 0;
                for (; p <= COL / mat_L; p += simd_inc)
                {
                    ele_vec = mat_simd_t::load(&ele[j][p], SIMD_ALIGN);
                    row_vec = mat_simd_t::load(&row[i - i_offset][p], SIMD_ALIGN);
                    row_vec ^= ele_vec;
                    xsimd::store(&row[i - i_offset][p], row_vec, SIMD_ALIGN);
                }
                for (; p <= COL / mat_L; p++)
                    row[i - i_offset][p] ^= ele[j][p];
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Send(row, n_workload * (LEN_LINE) * sizeof(mat_t),
             MPI_BYTE, 0, 1, MPI_COMM_WORLD);
}

int main()
{
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0)
    {
        mat_t(*ele)[LEN_LINE] =
            (mat_t(*)[LEN_LINE])aligned_alloc(128, sizeof(mat_t) * COL * (LEN_LINE));
        mat_t(*row)[LEN_LINE] =
            (mat_t(*)[LEN_LINE])aligned_alloc(128, sizeof(mat_t) * ROW * (LEN_LINE));

        memset(ele, 0, sizeof(mat_t) * COL * (LEN_LINE));
        memset(row, 0, sizeof(mat_t) * ROW * (LEN_LINE));

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

        timespec start, end;
        double time_used = 0;
        clock_gettime(CLOCK_REALTIME, &start);
        run_master(ele, row);
        clock_gettime(CLOCK_REALTIME, &end);
        time_used += (end.tv_sec - start.tv_sec) * 1000;
        time_used += double(end.tv_nsec - start.tv_nsec) / 1000000;
        cout << time_used << ',';

        free(ele);
        free(row);
    }
    else
    {
        run_slave();
    }

    MPI_Finalize();
}
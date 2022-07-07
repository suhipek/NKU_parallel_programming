#include <mpi.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <string.h>
#include <omp.h>

using namespace std;

#define ele_t float
#define ZERO (float)1e-5

#ifndef N
#define N 16
#endif

#ifndef DATA_PATH
#define DATA_PATH "/home/suhipek/NKU_parallel_programming/gauss.dat"
#endif

#define DEBUG

int world_size, world_rank;

void run_master(ele_t *_mat)
{
    int NUM_THREADS = world_size - 1;
    ele_t(*mat)[N] = (ele_t(*)[N])_mat;

    for (int i = 0; i < N; i++)
    {
        int n_lines = 0;
        if (world_size != 1)
        {
            n_lines = (N - i - 1) / (world_size - 1);
        }
        MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&n_lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(mat[i], sizeof(ele_t) * N, MPI_BYTE, 0, MPI_COMM_WORLD);
        for (int th = 1; th < world_size; th++)
        {
            MPI_Send(mat[i + 1 + (th - 1) * n_lines], sizeof(ele_t) * N * n_lines, MPI_BYTE, th, 0, MPI_COMM_WORLD);
        }
        int offset = NUM_THREADS ? NUM_THREADS * ((N - i - 1) / NUM_THREADS) : 0;
        for (int j = i + 1 + offset; j < N; j++)
        {
            if (mat[i][i] == 0)
                continue;
            ele_t div = mat[j][i] / mat[i][i];
            for (int k = i; k < N; k++)
                mat[j][k] -= mat[i][k] * div;
        }
        for (int th = 1; th < world_size; th++)
        {
            MPI_Recv(mat[i + 1 + (th - 1) * n_lines], sizeof(ele_t) * N * n_lines, MPI_BYTE, th, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    int i = -1;
    MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef DEBUG
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << mat[i][j] << ' ';
        cout << endl;
    }
    cout << endl;
#endif
}

void run_slave()
{
    int i_pivot;      // 枢轴位置
    int n_lines;      // 需要消元的行数
    ele_t lines_i[N]; // 当前消元行

    while (true)
    {
        MPI_Bcast(&i_pivot, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (i_pivot == -1)
            break;
        MPI_Bcast(&n_lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(lines_i, sizeof(ele_t) * N, MPI_BYTE, 0, MPI_COMM_WORLD);
        // printf("%d: pivot: %d, lines: %d\n", world_rank, i_pivot, n_lines);

        ele_t(*mat)[N] = (ele_t(*)[N])malloc(n_lines * N * sizeof(ele_t));
        MPI_Recv(mat, n_lines * N * sizeof(ele_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < n_lines; j++)
        {
            if (abs(mat[j][i_pivot]) < ZERO) // 枢轴为0，不需要消元
                continue;
            ele_t div = mat[j][i_pivot] / lines_i[i_pivot];
            for (int k = i_pivot; k < N; k++)
                mat[j][k] -= lines_i[k] * div;
        }
        free(mat);
        MPI_Send(mat, n_lines * N * sizeof(ele_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }
}

int main()
{
    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0)
    {
        ele_t *_mat = new ele_t[N * N];
        // ele_t(*mat)[N] = (ele_t(*)[N])_mat;
        ifstream data((string)DATA_PATH, ios::in | ios::binary);
        data.read((char *)_mat, N * N * sizeof(ele_t));
        data.close();

        timespec start, end;
        double time_used = 0;
        clock_gettime(CLOCK_REALTIME, &start);
        run_master(_mat);
        clock_gettime(CLOCK_REALTIME, &end);
        time_used += (end.tv_sec - start.tv_sec) * 1000;
        time_used += double(end.tv_nsec - start.tv_nsec) / 1000000;
        cout << time_used << endl;
        delete[] _mat;
    }
    else
    {
        run_slave();
    }

    MPI_Finalize();
}
#include <iostream>
#include <fstream>

#define N 1024

using namespace std;

int main()
{
    ofstream data("mar_vec.dat", ios::out | ios::binary);
    srand(314159265);
    int r;
    for(int i=0; i<(N*N+N); i++)
    {
        r = rand();
        data.write((char*)&r, sizeof(r));
    }
    data.close();
}
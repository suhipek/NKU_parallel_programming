#include <iostream>
#include <fstream>

#define N 128

using namespace std;

int main()
{
    ofstream data("gauss.dat", ios::out | ios::binary);
    srand(314159265);
    float r = 0;
    for (int i = 0; i < (N * N); i++)
    {
#ifdef R
        r = R;
#elif RP
        r = r + RP;
#else
        r = (float)rand();
#endif
        cout << r << ' ';
        data.write((char *)&r, sizeof(r));
    }
    data.close();
}
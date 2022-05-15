#include <iostream>
#include <fstream>

#define N 4096

using namespace std;

int main()
{
        ofstream data("gauss.dat", ios::out | ios::binary);
        srand(314159265);
        float r = 1;
#ifdef RC
        int rc = RC;
#endif
        for (int i = 0; i < (N * N); i++)
        {
#ifdef R
                r = R;
#elif RP
                r = r + RP;
#elif RC
                r = r + rc;
                rc++;
#else
                r = (float)rand();
#endif
                // cout << r << ' ';
                data.write((char *)&r, sizeof(r));
        }
        data.close();
}
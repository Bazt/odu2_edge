
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

//http://window.edu.ru/resource/958/40958/files/dvgu079.pdf
//http://mathhelpplanet.com/static.php?p=chislennyye-metody-resheniya-krayevykh-zadach
//http://blog.kislenko.net/show.php?id=1224
//https://www.cyberforum.ru/cpp-beginners/thread1160159.html

class solution {
public:
    double a, b, h;
    size_t n;
    
    vector<double> x;
    vector<double> p;
    vector<double> q;
    vector<double> f;
    
    double P(double x) {
        return 4*x;
    }

    double Q(double x) {
        return 4*x*x + 2;
    }

    double F(double x) {
        return 0;
    }
    
//    vector<vector<double>> A;
//    vector<vector<double>> B;

    
    vector<double> A, B, C, FF, k, v;
    
    void setup() {
        h = (b - a) / n;
        for (int i = 0; i <= n; i++) {
            x.push_back(a + i * h);
            p[i] = P(x[i]);
            q[i] = Q(x[i]);
            f[i] = F(x[i]);
        }
        
        A.resize(n - 1);
        B.resize(n - 1);
        C.resize(n - 1);
        FF.resize(n - 1);
        
        for (int i = 1; i <= n - 1; i++) {
            A[i] = 1 / h / h - p[i] / 2 / h;
            B[i] = 2 / h / h + q[i];
            C[i] = 1 / h / h + p[i] / 2 / h;
        }
        
        
        k.resize(n - 1);
        v.resize(n - 1);
        
        for (int j = 1; j <= n - 1; j++) {
            k[j] = B[j] / (C[j] + A[j] * k[j-1])
        }
        
    }
    
};

int main(int argc, const char * argv[]) {


    
    return 0;
}

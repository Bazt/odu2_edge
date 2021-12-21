
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

//http://window.edu.ru/resource/958/40958/files/dvgu079.pdf
//http://mathhelpplanet.com/static.php?p=chislennyye-metody-resheniya-krayevykh-zadach
//http://blog.kislenko.net/show.php?id=1224
//https://www.cyberforum.ru/cpp-beginners/thread1160159.html

class solution {
public:
    double a, b, h;
    double y0, yn;
    size_t n;
    
    vector<double> x;
    vector<double> p;
    vector<double> q;
    vector<double> f;
    
    vector<double> yy;
    
    double u_touch(double x) {
        return (x+1) * exp(-x*x);
    }
    
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
        x.resize(n+1);
        p.resize(n+1);
        q.resize(n+1);
        f.resize(n+1);
        for (int i = 0; i <= n; i++) {
            x[i] = (a + i * h);
            p[i] = P(x[i]);
            q[i] = Q(x[i]);
            f[i] = F(x[i]);
        }
        
        A.resize(n);
        B.resize(n);
        C.resize(n);
        FF.resize(n);
        
        for (int i = 1; i <= n - 1; i++) {
            A[i] = 1 / h / h - p[i] / 2 / h;
            C[i] = 2 / h / h - q[i];
            B[i] = 1 / h / h + p[i] / 2 / h;
            FF[i] = 0;
        }
        
        
        cout << "Check |C| > |A| + |B|: ";
        bool ok = true;
        for (int i = 1; i <= n - 1; i++) {
            if (abs(C[i]) > abs(A[i]) + abs(B[i])) {
                ok = false;
                break;
            }
        }
        cout << ok << endl;
        
        k.resize(n+1);
        v.resize(n+1);
        
        
        // b)
        k[0] = 0;
        v[0] = y0;
        
        k[n] = 0;
        v[n] = yn;
        
        for (int j = 1; j <= n - 1; j++) {
            k[j] = -B[j] / (C[j] + A[j] * k[j-1]);
            v[j] = - (A[j] * v[j - 1] - FF[j]) / (C[j] + A[j] * k[j-1]);
        }
        
        yy.resize(n+1);
        yy[n] = (v[n] + k[n]*v[n-1]) / (1 - k[n]*k[n-1]);
        for (int j = n-1; j >=0; j--) {
            yy[j] = k[j] * yy[j+1] + v[j];
        }
        
        
        
    }
    
};

int main(int argc, const char * argv[]) {
    solution s_b;
    s_b.a = 1e-3;
    s_b.b = 3;
    s_b.n = 100;
    s_b.y0 = s_b.u_touch(s_b.a);
    s_b.yn = s_b.u_touch(s_b.b);
    s_b.setup();
    for (int i = 0; i <= 100; i++) {
        cout << s_b.x[i] << "\t\t" << s_b.yy[i] << endl;
    }

    
    return 0;
}

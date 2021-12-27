
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

double u_tochnoe(double x) {
    return (x+1) * exp(-x*x);
}

double u_tochnoe_dx(double x) {
    return (1-2*x*x-2*x) * exp(-x*x);
}

class solution {
public:
    double a, b, h;
    double t_a, t_b;
    size_t n;
    double k1;
    double k2;
    double l1;
    double l2;
    
    solution(double a, double b, double k1, double k2, double l1, double l2, double t, double z, size_t n)
    : a(a), b(b), k1(k1), k2(k2), l1(l1), l2(l2), t_a(t), t_b(z), n(n) {    }
    
    
    vector<double> x;
    vector<double> p;
    vector<double> q;
    vector<double> f;
    
    vector<double> yy;

    
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
    
    
    void solve() {
        h = (b - a) / n;
        x.resize(n+1);
        p.resize(n+1);
        q.resize(n+1);
        f.resize(n+1);
        
        A.resize(n+1);
        B.resize(n+1);
        C.resize(n+1);
        FF.resize(n+1);
        
        for (int i = 0; i <= n; i++) {
            x[i] = a + i * h;
            A[i] = 1-P(x[i])*h/2;
            B[i] = 1 + P(x[i])*h/2;
            C[i] = 2-Q(x[i])*h*h;
            FF[i] = h*h*F(x[i]);
        }
        
        k.resize(n+1);
        v.resize(n+1);
        
        

        
        
        
        k[1] =-k2 / h*(k1-k2/h);
        v[1] = t_a / (k1 - k2/h);
        for (int j = 1; j <= n - 1; j++) {
            k[j+1] = B[j] / (C[j] - A[j] * k[j]);
            v[j+1] = (A[j] * v[j] - FF[j]) / (C[j] - A[j] * k[j]);
        }
        
        yy.resize(n+1);
        yy[n] = (l2*v[n] + t_b*h) / (l2+h*l1-l2*k[n]);
        for (size_t j = n; j >= 1; j--) {
            yy[j-1] = k[j] * yy[j] + v[j];
        }
        
        
        
    }
    
};

int main(int argc, const char * argv[]) {
    //solution(double a, double b, double k1, double k2, double l1, double l2, double t, double z, size_t n)
    
    // b) (4)
    ofstream f1("b4.txt");
    double a = 1e-3;
    double b = 3;
    double w_a = u_tochnoe(a);
    double w_b = u_tochnoe(b);
    solution s_b(a, b, 1, 0, 1, 0, w_a, w_b, 1000);
    s_b.solve();
    for (int i = 0; i <= s_b.n; i++) {
        // x_i << y_i
        f1 << s_b.x[i] << "\t\t" << s_b.yy[i] << endl;
    }

    // a) sposob 1
    
    ofstream f2("a1_1000.txt");
     a = 1e-3;
     b = 2;
     
    size_t n = 1000;
    double h = (b - a) / n;
    w_a = u_tochnoe_dx(a) * h;
    w_b = u_tochnoe(b);
    solution s_a1(a, b, -1, 1, 1, 0, w_a, w_b, n);
    s_a1.solve();
    for (int i = 0; i <= s_a1.n; i++) {
        // x_i << y_i
        f2 << s_a1.x[i] << "\t\t" << s_a1.yy[i] << endl;
    }
    
    ofstream f3("a1_10000.txt");
     a = 1e-3;
     b = 2;
     
     n = 10000;
     h = (b - a) / n;
    w_a = u_tochnoe_dx(a) * h;
    w_b = u_tochnoe(b);
    solution s_a1_1(a, b, -1, 1, 1, 0, w_a, w_b, n);
    s_a1_1.solve();
    for (int i = 0; i <= s_a1_1.n; i++) {
        // x_i << y_i
        f3 << s_a1_1.x[i] << "\t\t" << s_a1_1.yy[i] << endl;
    }
    return 0;
}

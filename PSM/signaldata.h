#ifndef SIGNALDATA_H
#define SIGNALDATA_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include<utils.h>

using namespace std;

template<typename T>
void print1d( T* in, int N){
    cout << fixed << setprecision(4);
    for(int ix = 0; ix < N; ix++){
        cout << in[ix] << " ";
    }
    cout << endl;
}

template<typename T>
void print2d( T* in, int N[2], int prec = 4, int cols = 0 ){
    if(cols !=  0){
        N[1] = cols;
    }
    cout << fixed << setprecision(prec);
    for(int iy = 0; iy < N[1]; iy++){
        cout << iy << ". ";
        for(int ix = 0; ix < N[0]; ix++){
            cout << in[iy + ix*N[1]] << " ";
        }
        cout << endl;
    }
    cout << endl;
}



template <typename T>
class point {
public:
    point(T a, T b, T c){
        x = a;
        y = b;
        z = c;
    }
    T absdiff(T a, T b, T c){
        return sqrt((x - a)*(x - a) + (y - b)*(y - b) + (z - c)*(z - c));
    }
    void showTargets(){
            cout << "x = " << x << " y = " << y  << " z = " << z <<endl;
    }
private:
    T x = 0;
    T y = 0;
    T z = 0;

};

class signaldata
{
public:
    signaldata(int n[3]);
    void addTarget(double x, double y, double z);
    void addTarget(point<double> p);
    void calcData(point <double> target);
    void showTargets();
    complex<double> *signalData;
    ~signaldata(){
        delete signalData;
        delete x;
    }
protected:
    double vel0 = 1.5*pow(10.0,8.0);
    double f = 8.95*pow(10.0,9.0);
    float *x;
    int N[3];
    vector<point<double>> targets;
};

#endif // SIGNALDATA_H

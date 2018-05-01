#include "signaldata.h"

signaldata::signaldata(int n[3])
{
    N[0] = n[0]; N[1] = n[1]; N[2] = n[2];
    x = new float[ N[1] ];
    x[0] = 0;
    for(int h = 1; h < N[1]; h++){
        x[h] = x[h-1] + 1.0/double(N[0]-1);
    }
    signalData = new complex<double>[N[0]*N[1]];
    for(int h = 0; h < N[0]*N[1]; h++){
        signalData[h] = 0;
    }
}

void signaldata::showTargets(){
    for(unsigned int h = 0; h < targets.size(); h++){
       targets[h].showTargets();
    }
}

void signaldata::addTarget(double x, double y, double z){
    point <double> p(x, y, z);
    targets.push_back(p);
    calcData(p);
}

void signaldata::addTarget(point<double> p){
    targets.push_back(p);
    calcData(p);
}

void signaldata::calcData(point <double> target){
    int timeIndizes[N[0]];
    for( int h = 0; h < N[0]; h++)
    {
       timeIndizes[h] = int(round(target.absdiff(x[h], 0.0, 0.0)*f/vel0));
    }
    for(int ix = 0; ix < N[0]; ix++){
        for( int it = 0; it < N[1]; it++){
            if(it == timeIndizes[ix] ){
                signalData[it + ix*N[1]] = signalData[it + ix*N[1]] + 1.0;
            }
        }
    }
}

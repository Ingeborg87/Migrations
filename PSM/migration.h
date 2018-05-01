#ifndef MIGRATION_H
#define MIGRATION_H
#include <signaldata.h>
#include <fftw3.h>

class migration:public signaldata
{
public:
    migration(int n[3]);
    void fillData(complex<double> *in);
    void showSizes();
    void calculateFac(double v);
    static void fftshift( int N[2], complex<double> *in, complex<double> *out, int dims, int dir);
    static void ifftshift( int N[2], complex<double> *in, complex<double> *out, int dims, int dir);
    static void multComplex(complex<double> *in, complex<double> *fac, int N );
    static void sum(complex<double> *in, complex<double> *out, int* N, int dim, int dir);
    void fftRawdata();
    void sumOverF();
    void ifftSummedData();
    void applyFac();
    void showRawData();
    void showImageData();
    void addZStepToImage(int nz);
    ~migration(){
        destroyAllWindows();
        delete kw; delete ekz; delete outx; delete kx; delete w; delete xz;
        fftw_destroy_plan(pfw);
        fftw_destroy_plan(pbw);
    }
private:
    int show = 0;
    double pi = acos(-1.0);
    double omega = 2.0*pi*f/vel0;
    double dz = 1.0/double(N[2]-1);
    double *kx;
    double *w;

    fftw_plan pfw;
    fftw_plan pbw;

    complex<double> *kw;
    complex<double> *ekz;
    complex<double> *xz;
    complex<double> *outx;


};

#endif // MIGRATION_H

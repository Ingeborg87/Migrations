#ifndef PSMUTILS_H
#define PSMUTILS_H
#include <iostream>
#include <algorithm>
#include <math.h>
#include <cstdint>
#include <cstring>
#include <fftw3.h>
template <class var>

class point
{
public:
    point(var x, var y, var z){
        this->x=x;
        this->y=y;
        this->z=z;
    }
    double absdiff(var x, var y, var z){
        return 2*sqrt(pow(this->x - x , 2) + pow(this->y - y, 2) + pow(this->z - z, 2));
    }

private:
    var x = 0.0;
    var y = 0.0;
    var z = 0.0;

};

void circshift(fftw_complex *in, fftw_complex*out, int xdim, int ydim, int xshift, int yshift);

void circshift1D_OP(fftw_complex *in, fftw_complex*out, int ydim, int yshift);

void ifftshift1D(fftw_complex *in, fftw_complex*out, int ydim);

void fftshift1D(fftw_complex *in, fftw_complex*out, int ydim);

void ifftshift2D(fftw_complex *in, fftw_complex*out, int xdim, int ydim);

void fftshift2D(fftw_complex *in, fftw_complex*out, int xdim, int ydim);


void circshift(fftw_complex *in, fftw_complex *out, int xdim, int ydim, int xshift, int yshift)
{
    if (xshift == 0 && yshift == 0)
    {
        out = in; //-- no change
        return;
    }

    for (int i = 0; i < xdim; i++)
    {
        int ii = (i + xshift) % xdim;
        if (ii < 0)
            ii = xdim + ii;
        for (int j = 0; j < ydim; j++)
        {
            int jj = (j + yshift) % ydim;
            if (jj < 0)
                jj = ydim + jj;
            out[ii * ydim + jj][0] = in[i * ydim + j][0];
                    out[ii * ydim + jj][1] = in[i * ydim + j][1];
        }
    }
}
 void circshift1D_IP(fftw_complex *in, int ydim, int yshift)
{
    if (yshift == 0)
        return;

    if (yshift > 0) // shift right
    {
        //std::rotate(&in[0], &in[ydim - yshift - 1], &in[ydim - 1]);
        std::rotate(in, in + (ydim - yshift), in + ydim);
    }
    else if (yshift < 0) // shift left
    {
        yshift = abs(yshift);
        //std::rotate(&in[0], &in[yshift], &in[ydim - 1]);
        std::rotate(in, in + yshift, in + ydim);
    }

    return;
}

void circshift1D_OP(fftw_complex *in, fftw_complex*out, int ydim, int yshift)
{
    if (yshift == 0)
    {
        out = in; //-- no change
        return;
    }

    memcpy(out, in, ydim * sizeof(fftw_complex));

    if (yshift > 0) // shift right
    {
        //std::rotate(&out[0], &out[ydim - yshift], &out[ydim]); // TODO check indices may be ydim-yshift
        std::rotate(out, out + (ydim - yshift), out + ydim); // C++ idiom: out + ydim is not used, out + ydim -1 is referenced
    }
    else if (yshift < 0) // shift left
    {
        yshift = abs(yshift);
        //std::rotate(&out[0], &out[yshift], &out[ydim - 1]);
        std::rotate(out, out + yshift, out + ydim); // TODO check
    }

    return;

//    for (int j = 0; j < ydim; j++)
//    {
//        int jj = (j + yshift) % ydim;
//        if (jj < 0)
//            jj = ydim + jj;
//        out[jj] = in[j];
//    }
}

void ifftshift1D(fftw_complex *in, fftw_complex*out, int ydim)
{
    //-- (ydim & 1)==0
    int pivot = (ydim % 2 == 0) ? (ydim / 2) : ((ydim + 1) / 2);
    //circshift1D(in, out, ydim, shiftBy);

    int rightHalf = ydim-pivot;
    int leftHalf = pivot;
    memcpy(out, in+(pivot), sizeof(fftw_complex)*rightHalf);
    memcpy(out+rightHalf, in, sizeof(fftw_complex)*leftHalf);
}

void fftshift1D(fftw_complex *in, fftw_complex*out, int ydim)
{
    int pivot = (ydim % 2 == 0) ? (ydim / 2) : ((ydim - 1) / 2);
    //circshift1D(in, out, ydim, shiftBy);
    int rightHalf = ydim-pivot;
    int leftHalf = pivot;
    memcpy(out, in+(pivot), sizeof(fftw_complex)*rightHalf);
    memcpy(out+rightHalf, in, sizeof(fftw_complex)*leftHalf);
}

void ifftshift2D(fftw_complex *in, fftw_complex*out, int xdim, int ydim)
{
    int shiftYBy = (ydim % 2 == 0) ? (ydim / 2) : ((ydim + 1) / 2);
    int shiftXBy = (xdim % 2 == 0) ? (xdim / 2) : ((xdim + 1) / 2);
    circshift(in, out, xdim, ydim, shiftXBy, shiftYBy);
}

void fftshift2D(fftw_complex *in, fftw_complex*out, int xdim, int ydim)
{
    int shiftYBy = (ydim % 2 == 0) ? (ydim / 2) : ((ydim - 1) / 2);
    int shiftXBy = (xdim % 2 == 0) ? (xdim / 2) : ((xdim - 1) / 2);
    circshift(in, out, xdim, ydim, shiftXBy, shiftYBy);
}

void multcomplex(fftw_complex *in, double phi){
    double phiIn = atan(in[0][0]/in[0][1]);
    double absIn = sqrt(in[0][0]*in[0][0] + in[0][1]*in[0][1]);
    in[0][0] = absIn*cos(phiIn + phi);
    in[0][1] = absIn*sin(phiIn + phi);
}

#endif // PSMUTILS_H

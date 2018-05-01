#include "migration.h"


migration::migration(int n[3]):signaldata(n)
{
    omega = 2.0*pi*f/vel0;
    dz = 1.0/double(N[2]-1);
    kx = new double[N[0]];
    w = new double[N[1]];
    kw = new complex<double>[N[0]*N[1]];
    ekz = new complex<double>[N[0]*N[1]] ;
    xz = new complex<double>[N[0]*N[2]];
    outx = new complex<double>[N[0]];
    pfw = fftw_plan_dft_2d(N[0],N[1],reinterpret_cast<fftw_complex*>(signalData),reinterpret_cast<fftw_complex*>(kw),FFTW_FORWARD,FFTW_ESTIMATE);
    pbw = fftw_plan_dft_1d(N[0],reinterpret_cast<fftw_complex*>(outx),reinterpret_cast<fftw_complex*>(outx),FFTW_BACKWARD,FFTW_ESTIMATE);
    if(N[0] % 2){
        kx[0] = -(N[0]-1)*pi;
    }
    else{
        kx[0] = -N[0]*pi;
    }
    if(N[1]%2){
        w[0]  = -double(N[1]-1)/double(2*(N[1] - 1));
    }
    else{
        w[0]  = -double(N[1])/double(2*N[1]);
    }
    for(int h = 1; h < N[0]; h++)
    {
        kx[h] = kx[h-1] + 2*pi;
    }
    for(int h = 1; h < N[1]; h++)
    {
        w[h] = w[h-1] + 1.0/double(N[1]-1.0);
    }
    calculateFac(vel0);
}

void migration::calculateFac(double v){
    vel0 = v;
    omega = 2.0*pi*f/vel0;
    double kz2 = 0.0;
    for(int it = 0; it < N[1]; it++){
        for( int ix = 0; ix < N[0]; ix++){
            kz2 = omega*omega*w[it]*w[it] - kx[ix]*kx[ix];
            if(kz2 >= 0.0){
                ekz[it + ix*N[1]] = exp(1i*sqrt(kz2)*dz);
            }
            else {
                ekz[it + ix*N[1]] = 1.0;
            }
        }
    }
}

void migration::showSizes(){
    cout << "t: "  << N[0] << endl;
    cout << "x: "  << N[1] << endl;
    cout << "z: "  << N[2] << endl;
}

void migration::fftshift( int N[2], complex<double> *in, complex<double> *out, int dims = 1, int dir = 0){
    // nsum - number of all elements in an 1d array
    if( sizeof(dir) != sizeof(int)){
        throw std::invalid_argument( "dir has to be a single integer!" );
    }
    if( dir < 0){
        throw std::invalid_argument( "dir has to be greater or equal than 0!" );
    }
    if(dir > 2){
        throw std::invalid_argument( "dir has to be smaller or equal than 2!" );
    }
    if( sizeof(dims) != sizeof(int)){
        throw std::invalid_argument( "dims has to be a single integer!" );
    }
    if( dims < 1){
        throw std::invalid_argument( "dims has to be greater or equal than 1!" );
    }
    if(dims > 2){
        throw std::invalid_argument( "dims has to be smaller or equal than 2!" );
    }
    if(dir > dims){
        throw std::invalid_argument( "dir has to be smaller or equal than dims!" );
    }
    int nSum = 1;
    for(int h = 0; h < dims; h++){
        nSum = nSum*N[h];
    }
    // create swap if in = out
    complex<double> *swap = new complex<double>[nSum];
    for(int h = 0; h < nSum; h++){
        swap[h] = in[h];
    }

    int nx, ny;
    if(dims==1){
        nx = N[0]; ny = 1;
    }
    else if(dims==2){
        nx = N[0]; ny = N[1];
    }
    int nShift = 0;
    if(dir == 1 || dir == 0){
        if(nx%2){
            nShift = (nx-1)/2;
            for(int hy = 0; hy < ny; hy++){
                for(int hx = nShift + 1; hx < nx; hx++){
                    out[hy + ny*(hx - nShift - 1 )] = swap[hy + ny*hx];
                    out[hy + ny*(hx - 1)] = swap[hy + ny*(hx - nShift - 1)];
                }
                out[hy + ny*(nx - 1)] = swap[hy + ny*nShift];
            }
        }
        else{
            nShift = nx/2;
            for(int hx = nShift; hx < nx; hx++){
                for(int hy = 0; hy < ny; hy++){
                    out[hy + ny*(hx - nShift)] = swap[hy + ny*hx];
                    out[hy + ny*hx] = swap[hy + ny*(hx - nShift)];
                }
            }
        }
    }
    if(dir == 2 || dir == 0){
        if(dir == 0){
            for(int h = 0; h < nSum; h++){
                swap[h] = out[h];
            }
        }
        if(ny%2){
            nShift = (ny-1)/2;
            for(int hx = 0; hx < nx; hx++){
                for(int hy = nShift + 1; hy < ny; hy++){
                    out[hy - nShift - 1 + ny*hx] = swap[hy + ny*hx];
                    out[hy -1 + ny*hx] = swap[hy - nShift - 1 + ny*hx];
                }
                out[ny - 1 + ny*hx] = swap[nShift + ny*hx];
            }
        }
        else{
            nShift = ny/2;
            for(int hy = nShift; hy < ny; hy++){
                for(int hx = 0; hx < nx; hx++){
                    out[hy - nShift + ny*hx] = swap[hy + ny*hx];
                    out[hy + ny*hx] = swap[hy - nShift + ny*hx];
                }
            }
        }
    }
    delete swap;
}

void migration::ifftshift( int N[2], complex<double> *in, complex<double> *out, int dims = 1, int dir = 0){
    // nsum - number of all elements in an 1d array
    if( sizeof(dir) != sizeof(int)){
        throw std::invalid_argument( "dir has to be a single integer!" );
    }
    if( dir < 0){
        throw std::invalid_argument( "dir has to be greater or equal than 0!" );
    }
    if(dir > 2){
        throw std::invalid_argument( "dir has to be smaller or equal than 2!" );
    }
    if( sizeof(dims) != sizeof(int)){
        throw std::invalid_argument( "dims has to be a single integer!" );
    }
    if( dims < 1){
        throw std::invalid_argument( "dims has to be greater or equal than 1!" );
    }
    if(dims > 2){
        throw std::invalid_argument( "dims has to be smaller or equal than 2!" );
    }
    if(dir > dims){
        throw std::invalid_argument( "dir has to be smaller or equal than dims!" );
    }
    int nSum = 1;
    for(int h = 0; h < dims; h++){
        nSum = nSum*N[h];
    }
    // create swap if in = out
    complex<double> *swap = new complex<double>[nSum];
    for(int h = 0; h < nSum; h++){
        swap[h] = in[h];
    }

    int nx, ny;
    if(dims==1){
        nx = N[0]; ny = 1;
    }
    else if(dims==2){
        nx = N[0]; ny = N[1];
    }
    int nShift = 0;
    if(dir == 1 || dir == 0){
        if(nx%2){
            nShift = (nx-1)/2;
            for(int hy = 0; hy < ny; hy++){
                for(int hx = nShift; hx < nx-1; hx++){
                    out[hy + ny*(hx - nShift)] = swap[hy + ny*hx];
                    out[hy + ny*(hx + 1)] = swap[hy + ny*(hx - nShift)];
                }
                out[hy + ny*nShift] = swap[hy + ny*(nx - 1)];
            }
        }
        else{
            nShift = nx/2;
            for(int hx = nShift; hx < nx; hx++){
                for(int hy = 0; hy < ny; hy++){
                    out[hy + ny*(hx - nShift)] = swap[hy + ny*hx];
                    out[hy + ny*hx] = swap[hy + ny*(hx - nShift)];
                }
            }
        }
    }
    if(dir == 2 || dir == 0){
        if(dir == 0){
            for(int h = 0; h < nSum; h++){
                swap[h] = out[h];
            }
        }
        if(ny%2){
            nShift = (ny-1)/2;
            for(int hx = 0; hx < nx; hx++){
                for(int hy = nShift; hy < ny-1; hy++){
                    out[hy - nShift + ny*hx] = swap[hy + ny*hx];
                    out[hy + 1 + ny*hx] = swap[hy - nShift + ny*hx];
                }
                out[nShift + ny*hx] = swap[ny - 1 + ny*hx];
            }
        }
        else{
            nShift = ny/2;
            for(int hy = nShift; hy < ny; hy++){
                for(int hx = 0; hx < nx; hx++){
                    out[hy - nShift + ny*hx] = swap[hy + ny*hx];
                    out[hy + ny*hx] = swap[hy - nShift + ny*hx];
                }
            }
        }
    }
    delete swap;
}

void migration::multComplex(complex<double> *in, complex<double> *fac, int N ){
    for(int h = 0; h < N; h++){
        in[h] = in[h]*fac[h];
    }

}

void migration::sum(complex<double> *in, complex<double> *out, int* N, int dim, int dir){
    if(dim == 1){
        out[0] = 0.0;
        for(int h = 0; h < N[0]; h++){
            out[0] = out[0] + in[h];
        }
    }
    if(dim == 2){
        if(dir == 1){
            for(int h = 0; h < N[1]; h++){
                out[h] = 0.0;
            }
            for(int hy = 0; hy < N[1]; hy++){
                for(int hx = 0; hx < N[0]; hx++){
                    out[hy] = out[hy] + in[hy + N[1]*hx];
                }
            }
        }
        if(dir == 2){
            for(int h = 0; h < N[0]; h++){
                out[h] = 0.0;
            }
            for(int hx = 0; hx < N[0]; hx++){
                for(int hy = 0; hy < N[1]; hy++){
                    out[hx] = out[hx] + in[hy + N[1]*hx];
                }
            }
        }
    }
}

void migration::sumOverF(){
    int n[2] = {N[0], N[1]};
    sum(kw, outx, n,2,2);
}

void migration::fftRawdata(){
    int n[2];
    n[0] = N[0]; n[1] = N[1];
    fftw_execute(pfw);
    fftshift(n, kw, kw, 2,0);
}

void migration::ifftSummedData(){
    int n[2];
    n[0] = N[0]; n[1] = N[1];
    ifftshift(n, outx, outx, 1,0);
    fftw_execute(pbw);
}

void migration::applyFac(){
    multComplex(kw, ekz, N[0]*N[1]);
}

void migration::showRawData(){
    int Nz[2];
    Nz[0] = N[0];
    Nz[1] = N[2];
    Mat image = Mat::zeros(2,Nz,CV_8UC1), coloredImage;
    double maxVal = real(signalData[0]), minVal = real(signalData[0]);
    for(int h = 0; h < Nz[0]*Nz[1]; h++){
        if(maxVal < real(signalData[h]) ){
            maxVal = real(signalData[h]);
        }
        if(minVal > real(signalData[h])){
            minVal = real(signalData[h]);
        }
    }
    for ( int iz=0;iz<Nz[1];iz++) {
        for ( int ix=0;ix<Nz[0];ix++) {
            image.at<uchar>(iz,ix) =  uchar(round(255.0*(real(signalData[iz + ix*Nz[1]]) - minVal)/(maxVal-minVal)));
        }
    }
    namedWindow("rawdata", WINDOW_FREERATIO);
    applyColorMap(image, coloredImage, COLORMAP_JET);
    imshow("rawdata", coloredImage );
}

void migration::showImageData(){
    int Nz[2];
    Nz[0] = N[0];
    Nz[1] = N[2];
    Mat image = Mat::zeros(2,Nz,CV_8UC1), coloredImage;
    double maxVal = real(xz[0]), minVal = real(xz[0]);
    for(int h = 0; h < Nz[0]*Nz[1]; h++){
        if(maxVal < real(xz[h]) ){
            maxVal = real(xz[h]);
        }
        if(minVal > real(xz[h]) ){
            minVal = real(xz[h]);
        }
    }
    for ( int iz=0;iz<Nz[1];iz++) {
        for ( int ix=0;ix<Nz[0];ix++) {
            image.at<uchar>(iz,ix) =  uchar(round(255.0*(real(xz[iz + ix*Nz[1]]) - minVal)/(maxVal-minVal)));
        }
    }
    namedWindow("imagedata", WINDOW_FREERATIO);
    applyColorMap(image, coloredImage, COLORMAP_JET);
    imshow("imagedata", coloredImage );
}

void migration::addZStepToImage(int nz){
    for(int ix = 0; ix < N[0]; ix++){
        xz[nz + ix*N[2]] = outx[ix];
    }
}

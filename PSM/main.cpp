
#include <signaldata.h>
#include <opencv2/opencv.hpp>
#include <migration.h>
using namespace cv;

int main()
{
    /// setting up modeldata
    /// // velovities, frequences etc
    int sizesX = 501, sizesZ = 501, sizesT = 501;
    int N[3] = {sizesX, sizesT, sizesZ};
    point <double> p(0.15,0.0,0.75);
    migration psm(N);
    psm.addTarget(p);
    psm.addTarget(0.7,0,0.2);
    psm.addTarget(0.5,0,0.5);
    psm.addTarget(0.15,0,0.25);
    psm.addTarget(0.75, 0, 0.75);
    psm.showRawData();
    psm.fftRawdata();

    for(int iz = 0; iz < sizesZ; iz++){
        psm.applyFac();
        psm.sumOverF();
        psm.ifftSummedData();
        psm.addZStepToImage(iz);
        cout << "Shift: " << iz + 1 << "/" << N[2] << endl;
    }
    psm.showImageData();
    waitKey();
}


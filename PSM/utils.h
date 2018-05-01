#ifndef UTILS
#define UTILS
#include <iostream>
#include <iomanip>
#include <complex>
#include <opencv2/opencv.hpp>
#include <fftw3.h>

using namespace std;
using namespace cv;


template<typename T>
void showData( T* in, int sizes[2]){
    Mat image = Mat::zeros(2,sizes,CV_8UC1), coloredImage;
    double maxVal = in[0], minVal = in[0];
    for(int h = 0; h < sizes[0]*sizes[1]; h++){
        if(maxVal < in[h] ){
            maxVal = in[h];
        }
        if(minVal > in[h] ){
            minVal = in[h];
        }
    }
    for ( int iz=0;iz<sizes[1];iz++) {
        for ( int ix=0;ix<sizes[0];ix++) {
            image.at<uchar>(iz,ix) =  uchar(round(255.0*(in[iz + ix*sizes[1]] - minVal)/(maxVal-minVal)));
        }
    }
    namedWindow("Image", WINDOW_FREERATIO);
    applyColorMap(image, coloredImage, COLORMAP_HOT);
    imshow("Image", coloredImage );
    waitKey();
    destroyAllWindows();
}

#endif // UTILS


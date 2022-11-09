#pragma once
#include <complex>
#include <vector>
#include "typedf.h"

using namespace std;
class Encoder{
    public:
        int N;
        double scale;
        dcomplex xi;
        complexArray sigma_R_basis;
        complexMatrix FFTMatrix;
        complexMatrix IFFTMatrix;  
        Encoder(int N,double scale);
        complexMatrix vandermonde(dcomplex xi);
        complexMatrix vandermondeInv(dcomplex xi);
        complexArray sigma(complexArray x);
        complexArray sigmaInv(complexArray x); 
        IntegerArray encode(complexArray x);
        complexArray decode(IntegerArray x);
};

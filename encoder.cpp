#include "encoder.h"
#include <iostream>
#include <cmath>

Encoder::Encoder(int N, double scale) {
    this->N = N;
    this->scale = scale;
    this->xi = complex<double>(cos(2 * M_PI / (2 * N)),sin(2 * M_PI / (2 * N)));
    this->scale = scale;
    this->FFTMatrix = Encoder::vandermonde(xi);
    this->IFFTMatrix = Encoder::vandermondeInv(xi);
}

complexMatrix Encoder::Encoder::vandermonde(dcomplex xi) {
    cout << N << endl;
    complexMatrix matrix(N,complexArray(N));
    for(int i = 0; i < N; i++){
        dcomplex root = pow(xi,i);
        for(int j = 0; j < N; j++){
            matrix[j][i] = pow(root,(2 * j + 1));
        }
    }
    return matrix; 
}

complexMatrix Encoder::Encoder::vandermondeInv(dcomplex xi) {
    complexMatrix matrix(N,complexArray(N));
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            matrix[j][i] = this->FFTMatrix[i][j] / dcomplex(N);
        }
    }
    return matrix; 
}

complexArray Encoder::sigma(complexArray x){
    complexArray res(N);
    for(int i = 0; i < N; i++){
        dcomplex temp(0,0);
        for(int j = 0; j < N; j++){
            temp += FFTMatrix[i][j] * x[j];
        }
        res[i] = temp;
    }
    return res;
}

complexArray Encoder::sigmaInv(complexArray x){
    complexArray res(N);
    for(int i = 0; i < N; i++){
        dcomplex temp(0,0);
        for(int j = 0; j < N; j++){
            temp += IFFTMatrix[i][j] * x[j];
        }
        res[i] = temp;
    }
    return res;
}

IntegerArray Encoder::encode(complexArray x){
    if(x.size() != N/2){
        throw "size error";
    }
    for(int i = 0; i < N/2; i++){
        x[i] = x[i] * dcomplex(scale);
    }
    for(int i = N/2 - 1; i >= 0; i--){
        x.push_back(conj(x[i]));
    }

    complexArray temp = sigmaInv(x);
    // for(int i = 0; i < N; i++){
    //     cout << temp[i].real() << "+" << temp[i].imag() << "i\t";
    // }
    // cout << endl;

    IntegerArray res(N);
    for(int i = 0; i < N;i++){
        res[i] = (long long)(temp[i].real());
    }
    return res;
}
complexArray Encoder::decode(IntegerArray x){
    complexArray res(N);
    for(int i = 0; i < N; i++){
        res[i] = dcomplex(x[i]*1.0/scale);
    }
    res = sigma(res);
    complexArray ans;
    for(int i = N - 1;i >= N/2;i--){
        ans.push_back(res[i]);
    }
    return ans;
}
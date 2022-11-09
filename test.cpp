#include "encoder.h"
#include <iostream>
using namespace std;
const int N = 4;
int main(){
    
    Encoder encoder(N,1024);
    complexArray x(N/2);
    x[0] = dcomplex(1.23,1);
    x[1] = dcomplex(2.14,2);
    // x[2] = dcomplex(3.23,1);
    // x[3] = dcomplex(4.14,2);
    for(int i = 0; i < N/2; i++){
        cout << x[i].real() << "+" << x[i].imag() << "i\t";
    }
    cout << endl;

    // IntegerArray p = encoder.encode(x);
    IntegerArray p(N,1);
    
    for(int i = 0; i < N; i++){
        cout << p[i] << " ";
    }
    cout << endl;
    complexArray y = encoder.decode(p);
    for(int i = 0; i < N/2; i++){
        cout << y[i].real() << "+" << y[i].imag() << "i\t";
    }
    cout << endl;
}
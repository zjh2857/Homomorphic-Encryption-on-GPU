#include "./ckks.h"
#include <ctime>

ciphertext::ciphertext(IntegerArray plaintext) {
    this->s = randompoly(plaintext.size(),1000);
    this->e = randompoly(plaintext.size(),10);
    IntegerArray t = polymul(plaintext,s);
    this->c = polyadd(t,e);       
}

ciphertext ciphertext::mulplaintext(IntegerArray plaintext) {
    
}

IntegerArray ciphertext::polymul(IntegerArray a, IntegerArray b) {
    
}

IntegerArray ciphertext::polyadd(IntegerArray a, IntegerArray b) {
    
}

IntegerArray ciphertext::randompoly(int len,long long upper){
    IntegerArray res(len);
    srand((int)time(0));
    for(int i = 0; i < len; i++){
        res[i] = rand()%upper;
    }
    return res;
} 
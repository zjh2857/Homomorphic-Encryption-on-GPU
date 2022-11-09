#pragma once
#include "typedf.h"

class ciphertext{
    public:
        IntegerArray c;
        IntegerArray s;
        ciphertext(IntegerArray plaintext);
        ciphertext mulplaintext(IntegerArray plaintext);
    private:
        IntegerArray e;
        IntegerArray polymul(IntegerArray a,IntegerArray b);
        IntegerArray polyadd(IntegerArray a,IntegerArray b);
        IntegerArray randompoly(int len,long long upper);    
};

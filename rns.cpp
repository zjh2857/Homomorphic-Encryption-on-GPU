#define N 8
#include <cstdio>
#include <stdlib.h>
// #include "uint128_t.h"
// #include <bits/stdc++.h>
using namespace std;
unsigned long long qpow(unsigned long long a,unsigned long long b,unsigned long long q){
    unsigned long long r=1;
    // unsigned long long base = a;
    while(b){
        if(b&1)r = (a * r)%q;
        a = a * a % q;
        b >>= 1;
    }
    return r;
}
void print(unsigned long long* a){
    for(int i = 0; i < N; i++){
        printf("%llu\t",a[i]);
    }
    printf("\n");
}
void bigIntegerMul(unsigned long long * a, unsigned long long b,int size){
    __uint128_t temp = a[0];
    temp = temp * b;
    a[0] = temp;
    unsigned long long carry = temp >> 64;
    
    for(int i = 1; i < size; i++){
        if(!a[i] && !carry){
            break;
        }
        __uint128_t temp = a[i];
        temp = temp * b;
        a[i] = (temp & 0xffffffffffffffff);
        if(a[i] + carry < a[i]){
            a[i] += carry;
            carry = 1 + temp >> 64;
        }
        else{
            a[i] += carry;
            carry = 1 + temp >> 64;
        }
    }
}
void bigIntegerAdd(unsigned long long * a, unsigned long long *b,int size){
    unsigned long long carry = 0;
    for(int i = 0; i < size; i++){
        if(a[i] + b[i] + carry < a[i]){
            // carry = 1;
            a[i] = a[i] + b[i] + carry;
            carry = 1;
        }
        else{
            a[i] = a[i] + b[i] + carry;
            carry = 0;
        }
    }
}
int isneg(unsigned long long * a, unsigned long long *b ,unsigned long long p,int size){
    unsigned long long* temp = (unsigned long long*)malloc(size * sizeof(unsigned long long)); 
    for(int i = 0; i < size; i++){
        temp[i] = a[i];
    }
    bigIntegerMul(temp,p,size);
    for(int i = 0; i < size - 1; i++){
        if(a[i] > temp[i]){
            free(temp);
            return 1;
        }
        else if(a[i] < temp[i]){
            free(temp);
            return -1;
        }
    }
    free(temp);
    return 0;
}
unsigned long long bigIntegerMod(unsigned long long * a, unsigned long long *b,int size){
    unsigned long long l = 0;
    unsigned long long r = (1 << 63);
    while(l <= r){
        unsigned long long guess = (l + r)/2;
        int res = isneg(a,b,guess,size);
        if(res == 1){
            l = guess + 1;
        }
        else if(res == -1){
            r = guess;
        }
        else{
            return a[0] - (b[0] * guess);
        }
    }
    return 0;
}
int main(){
    unsigned long long moduleChain[N] = {10007,10009,10037,10039,10061,10067,10069,10079};
    unsigned long long Ni[N][N];
    unsigned long long bigN[N];
    for(int i = 0; i < N;i++){
        for(int j = 0; j < N; j++){
            Ni[i][j] = 0;
        }
    }
    for(int i = 0; i < N;i++){
        bigN[i] = 0;
    }
    bigN[0] = 1;
    unsigned long long ti[N];
    unsigned long long bignum1 = 11451;
    unsigned long long bignum1Rns[N];
    unsigned long long bignum2 = 37777;
    unsigned long long bignum2Rns[N];
    for(int i = 0; i < N; i++){
        bignum1Rns[i] = bignum1 % moduleChain[i];
        // printf("%llu,",bignum1Rns[i]);
        bignum2Rns[i] = bignum2 % moduleChain[i];
    }
    for(int i = 0;i < N; i++){
        Ni[i][0] = 1;
        ti[i] = 1;
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(i==j)continue;
            bigIntegerMul(Ni[i],moduleChain[j],N);
            ti[i] = ti[i] * qpow(moduleChain[j],moduleChain[i]-2,moduleChain[i]) % moduleChain[i];
        }
    }
    
    for(int i = 0; i < N; i++){
        bigIntegerMul(bigN,moduleChain[i],N);
    }
    unsigned long long recovery = 0;
    unsigned long long res[N];
    for(int i = 0; i < N; i++)res[i] = 0;
    for(int i = 0; i < N; i++){
        bigIntegerMul(Ni[i],ti[i],N);
        bigIntegerMul(Ni[i],bignum1Rns[i],N);
        bigIntegerAdd(res,Ni[i],N);
    }
    print(bigN);
    print(res);
    recovery = bigIntegerMod(res,bigN,N);
    
    printf("%llu\n",recovery);
}
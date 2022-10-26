#include <iostream>
using namespace std; 
long long M = 998244353;
long long g = 3;

long long qpow(long long a,long long b,long long n){
    long long ans = 1;
    long long base = a;
    while(b){
        if(b&1){
            ans *= base;
            ans %= n;
        }
        base = base * base % n;
        b >>= 1;
    }
    return ans;
}

int main(){
    cout << qpow(2,7,5);
}
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

const long long M = 998244353;
const long long N = 32;
long long Ninv;
long long revphi[N/2];
long long revphiinv[N/2];
long long bitrevs[N];
long long bitrevl[N];
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

long long _bitrevs(long long x){
    long long res = 0;
    int l = 0;
    long long temp = N;
    while(temp != 1){
        temp=temp>>1;
        l++;
    }
    for(int i = 0; i < l - 1; i++){
        // cout << x << endl;
        res = (res << 1) | (x & 1);
        x >>= 1;
    }
    return res;
}
long long _bitrevl(long long x){
    long long res = 0;
    int l = 0;
    long long temp = N;
    while(temp != 1){
        temp=temp>>1;
        l++;
    }
    for(int i = 0; i < l; i++){
        // cout << x << endl;
        res = (res << 1) | (x & 1);
        x >>= 1;
    }
    return res;
}
void init(){
    long long g = qpow(3,(M-1)/N,M);
    long long gi = qpow(332748118,(M-1)/N,M);
    for(int i = 0; i < N; i++){
        bitrevs[i] = _bitrevs(i);
        bitrevl[i] = _bitrevl(i);
    }
    for(int i = 0; i < N/2; i++){
        revphi[i] = qpow(g,bitrevs[i],M);
        revphiinv[i] = qpow(gi,bitrevs[i],M);
    }
    Ninv = qpow(N,M-2,M);
}
void swapInv(long long* a){
    for(int i = 0; i < N; i++){
        if(i < bitrevl[i])swap(a[i],a[bitrevl[i]]);
    }
    // for(int i = 0; i < N; i++){
    //     cout << a[i] << " ";
    // }
    // cout << endl;
}
long long* NTT(long long* a){
    long long t = N;
    long long m = 1;
    while(t > 1){
        t /= 2;
        for(int i = 0; i < m; i++){
            int j1 = 2 * i * t;
            int j2 = j1 + t;
            for(int j = j1; j < j2; j++){
                long long u = a[j] % M;
                long long v = a[j + t] * revphi[i] % M;
                a[j] = (u + v + M) % M;
                a[j + t] = (u - v + M) % M;
            }
        }
        m <<= 1;
    }
    return a;
}

long long* INTT(long long* a){
    swapInv(a);
    long long t = N;
    long long m = 1;
    while(t > 1){
        t /= 2;
        for(int i = 0; i < m; i++){
            int j1 = 2 * i * t;
            int j2 = j1 + t;
            for(int j = j1; j < j2; j++){
                long long u = a[j] % M;
                long long v = a[j + t] * revphiinv[i] % M;
                a[j] = (u + v + M) % M;
                a[j + t] = (u - v + M) % M;
            }
        }

        m <<= 1;
    }
    for(int i = 0; i < N; i++){
        a[i] = a[i] * Ninv % M;
    }
    swapInv(a);
    return a;
}
// __global__ void NTT(long long* a){
//     int i = threa
// } 

int main(){
    long long a[N];
    for(int i = 0; i < N;i++){
        a[i] = i*1235%100000;
    }
    for(int i = 0; i < N; i++){
        cout << a[i] << " ";
    }
    cout << endl;
    long long b[N];
    for(int i = 0; i < N;i++){
        b[i] = i*2;
    }
    long long c[N];
    init();
    INTT(a);
    // INTT(b);
    srand((int)time(0)); 
    for(int i = 0; i < N; i++){
        // cout << a[i] <<"," << b[i] << endl;
        c[i] = (a[i]+(rand()%2)-1)%M;
    }
    NTT(c);
    for(int i = 0; i < N; i++){
        cout << c[i] << " ";
    }
}
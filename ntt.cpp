#include<bits/stdc++.h>
#define poi 2100000
using namespace std;
typedef long long ll;
const int mod=998244353,yg=3;//模数 998244353,原根为3 
int a[poi],b[poi],lim=1,l,r[poi];
inline int re() {
	char x = getchar();
	int k = 1, y = 0;
	while(x < '0' || x > '9')
	{if(x == '-') k = -1;x = getchar();}
	while(x >= '0' && x <= '9')
	{y = (y << 3) + (y << 1) + x - '0'; x = getchar();}
	return y * k;
}
inline void wr(int x) {
	if(x < 0) putchar('-'), x = -x;
	if(x > 9) wr(x / 10);
	putchar(x % 10 + '0');
}
inline int ksm(int x,int y) {
	int ans = 1;
	for(; y; y >>= 1)
	{
		if(y & 1) ans = (ll)((ll)ans * x) % mod;
		x = (ll)((ll)x * x) % mod;
	}
	return ans;
}
void ntt(int *a, bool tp)
{
	for(int i = 0; i < lim; i++)
	if(i < r[i]) swap(a[i], a[r[i]]);
	for(int mid = 1; mid < lim; mid <<= 1)
	{
		int bas = ksm(tp ? yg : 332748118, (mod - 1) / (mid << 1));
		for(int i = mid << 1, j = 0; j < lim; j += i)
		{
			int w = 1;
			for(int k = 0; k < mid; k++, w = ((ll)w * bas) % mod)
			{
				int x = a[k + j], y = (ll)w * (ll)a[k + j + mid] % (ll)mod;
				a[k + j] = (x + y) % mod;
				a[k + mid + j] = (x - y + mod) % mod;
			}
		}
	}
}
signed main() {
	int n = 4, m = 4;
	for(int i = 0; i <= n; i++) a[i] = i;
	for(int i = 0; i <= m; i++) b[i] = 2 * i;
    // b[3] = 1;
    // a[2] = 1;
	while(lim < 4) lim <<= 1, l++;
	for(int i = 0; i < lim; i++)
	r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));
	ntt(a, 1), ntt(b, 1);
	for(int i = 0; i < lim; i++){
		printf("%llu,%llu\n",a[i],b[i]);
	}
	for(int i = 0; i < lim; i++) a[i] = ((ll)a[i] + b[i]) % mod;
	ntt(a, 0);
	int niyuan = ksm(lim, mod-2);
	for(int i = 0;i <= 4; i++)
	wr(((ll)a[i] * niyuan) % mod), putchar(' ');
	return 0;
}
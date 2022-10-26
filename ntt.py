import numpy as np

M = 998244353
g = 3
N = 8
g = pow(g,((M-1) // N),M)

phi = np.zeros((N,N),np.int64)
for i in range(N):
    for j in range(N):
        phi[i][j] = pow(g,i*j,M)

invphi = np.zeros((N,N),np.int64)
for i in range(N):
    for j in range(N):
        invphi[i][j] = pow(g,-i*j,M)

invphi = (invphi*pow(N,-1,M)%M)
revphi = np.zeros(N,np.int64)
for i in range(N//2):
    revphi[i] = pow(g,int(bin(i)[2:][::-1].ljust(len(bin(N))-4,"0"),2),M)
rev = [0] * N
for i in range(N):
    rev[i] = int(bin(i)[2:][::-1].ljust(len(bin(N))-3,"0"),2)
# print(rev)
def ntt(aa,inv):
    global g
    
    if(inv):
        g = pow(g,-1,M)
        b = np.zeros(N,np.int64)
        for i in range(len(aa)):
            b[i] = aa[rev[i]] 
        a = b

    else:
        a = aa
    t = N
    m = 1
    while(t > 1):
        t //= 2
        for i in range(m):
            j1 = 2 * i * t
            j2 = j1 + t
            for j in range(j1,j2):
                u = a[j]
                v = a[j+t] * revphi[i] % M
                a[j] = (u + v) % M
                a[j + t] = (u - v) % M
        m *= 2
    return np.array(a,np.int64)

# a = np.array([0,0,0,0,0,0,1,2,0,0,0,0,0,0,1,2]).T
# b = np.array([0,0,0,0,0,1,2,1,0,0,0,0,0,0,1,2]).T
# a = phi.dot(a) % M 
# b = phi.dot(b) % M
# t = a * b % M
# print(invphi.dot(t)%M)

a = np.array([0,0,0,0,0,0,1,2]).T
b = np.array([0,0,0,0,0,0,1,2]).T
a = ntt(a,False) % M
# print(revphi)
print(a)
b = ntt(b,False) % M
t = a * b % M
# print(ntt(t,True) * pow(N,-1,M) % M)


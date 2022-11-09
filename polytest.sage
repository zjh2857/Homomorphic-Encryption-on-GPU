R.<t> = PolynomialRing(Zmod(7)) 
m = t ^ 3 + t ^ 2 
a = 2 * t ^ 2 + 5
n = t ^ 4 + 1

s = m * a % n

def exgcd(a,b):
    if(b == 0*t):
        return 1,0,a
    else:
        x,y,q = exgcd(b,a%b)
        x,y = y, (x - (a // b) * y)
        return x,y,q
a = t ^ 5 
e = 1 + t
n = t^4 + 1
print(a % n)
# a_inv = exgcd(a,n)[0]

# print(a_inv)
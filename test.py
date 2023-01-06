import numpy as np
M = 32
xi = np.exp(2 * np.pi * 1j / M)
yi = np.exp(-2 * np.pi * 1j / M)
def vandermonde(xi: np.complex128, M: int) -> np.array:
    """Computes the Vandermonde matrix from a m-th root of unity."""
    
    N = M //2
    matrix = []
    # We will generate each row of the matrix
    for i in range(N):
        # For each row we select a different root
        root = xi ** (2 * i + 1)
        row = []

        # Then we store its powers
        for j in range(N):
            row.append(root ** j)
        matrix.append(row)
    return np.array(matrix)

def vandermondeinv(xi: np.complex128, M: int) -> np.array:
    """Computes the Vandermonde matrix from a m-th root of unity."""
    
    N = M //2
    matrix = []
    # We will generate each row of the matrix
    for i in range(N):
        # For each row we select a different root
        root = xi ** (2 * i + 1)
        row = []

        # Then we store its powers
        for j in range(N):
            row.append(root ** (-j))
        matrix.append(row)
    return np.array(matrix)
for i in range(9):
    print(np.round(pow(xi,i)*10000))
A = vandermonde(xi,M)
C = vandermondeinv(xi,M).T / (M//2)
z = [0] * (M//2)
z[0] = 100
z[-1] = 100 
z = np.array(z)
t = (C).dot(z)
print(t)
# res = A.T.dot(C)
# print(np.round(res*10000))
# A = vandermonde(xi,M)
# C = vandermondeinv(xi,M)
B = np.linalg.inv(A)

# # print(C)
# # print(B.T*8-A)
# # z = [0] * (M//2)
# # z[0] = 100
# # z[-1] = 100
# # z = np.array(z)
# # t = (A.T).dot(z)
# # print(t)
# z = np.array([1,1,3,4])
# t = np.array([1,0,0,0])
# zz = B.dot(z)
# tt = (A.dot(t))

# # print(A.)
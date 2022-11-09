import numpy as np
M = 8
xi = np.exp(2 * np.pi * 1j / M)

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

A = vandermonde(xi,M)

B = np.linalg.inv(A)
for i in range(9):
    print(i,pow(xi,i))
print(A)
print(B.T*4)
# z = np.array([3 +4j, 2 - 1j,2 + 1j,3-4j])
# t = B.dot(z)
# print(t)
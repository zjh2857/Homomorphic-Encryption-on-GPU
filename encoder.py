#引入包
# !pip3 install fastcore
import numpy as np
from numpy.polynomial import Polynomial
from fastcore.foundation import patch_to

#CKKSEncoder类
class CKKSEncoder:
    """Basic CKKS encoder to encode complex vectors into polynomials."""
    
    def __init__(self, M: int):
        """Initialization of the encoder for M a power of 2. 
        
        xi, which is an M-th root of unity will, be used as a basis for our computations.
        """
        self.xi = np.exp(2 * np.pi * 1j / M)
        self.M = M
        
    @staticmethod
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
        return matrix
    
    def sigma_inverse(self, b: np.array) -> Polynomial:
        """Encodes the vector b in a polynomial using an M-th root of unity."""

        # First we create the Vandermonde matrix
        A = CKKSEncoder.vandermonde(self.xi, M)

        # Then we solve the system
        coeffs = np.linalg.solve(A, b)

        # Finally we output the polynomial
        p = Polynomial(coeffs)
        return p

    def sigma(self, p: Polynomial) -> np.array:
        """Decodes a polynomial by applying it to the M-th roots of unity."""

        outputs = []
        N = self.M //2

        # We simply apply the polynomial on the roots
        for i in range(N):
            root = self.xi ** (2 * i + 1)
            output = p(root)
            outputs.append(output)
        return np.array(outputs)

#\pi 映射和 \pi逆映射
@patch_to(CKKSEncoder)
def pi(self, z: np.array) -> np.array:
    """Projects a vector of H into C^{N/2}."""
    
    N = self.M // 4
    return z[:N]

@patch_to(CKKSEncoder)
def pi_inverse(self, z: np.array) -> np.array:
    """Expands a vector of C^{N/2} by expanding it with its
    complex conjugate."""
    
    z_conjugate = z[::-1]
    z_conjugate = [np.conjugate(x) for x in z_conjugate]
    return np.concatenate([z, z_conjugate])

#生成本原根
@patch_to(CKKSEncoder)
def create_sigma_R_basis(self):
    """Creates the basis (sigma(1), sigma(X), ..., sigma(X** N-1))."""

    self.sigma_R_basis = np.array(self.vandermonde(self.xi, self.M)).T

#
@patch_to(CKKSEncoder)
def compute_basis_coordinates(self, z):
    """Computes the coordinates of a vector with respect to the orthogonal lattice basis."""
    output = np.array([np.real(np.vdot(z, b) / np.vdot(b,b)) for b in self.sigma_R_basis])
    return output

def round_coordinates(coordinates):
    """Gives the integral rest."""
    coordinates = coordinates - np.floor(coordinates)
    return coordinates

def coordinate_wise_random_rounding(coordinates):
    """Rounds coordinates randonmly."""
    r = round_coordinates(coordinates)
    f = np.array([np.random.choice([c, c-1], 1, p=[1-c, c]) for c in r]).reshape(-1)
    
    rounded_coordinates = coordinates - f
    rounded_coordinates = [int(coeff) for coeff in rounded_coordinates]
    return rounded_coordinates

@patch_to(CKKSEncoder)
def sigma_R_discretization(self, z):
    """Projects a vector on the lattice using coordinate wise random rounding."""
    coordinates = self.compute_basis_coordinates(z)
    print("coordinates",coordinates)
    rounded_coordinates = coordinate_wise_random_rounding(coordinates)
    print("rounded_coordinates",rounded_coordinates)
    y = np.matmul(self.sigma_R_basis.T, rounded_coordinates)
    return y

#编码和解码
@patch_to(CKKSEncoder)
def __init__(self, M:int, scale:float):
    """Initializes with scale."""
    self.xi = np.exp(2 * np.pi * 1j / M)
    self.M = M
    self.create_sigma_R_basis()
    self.scale = scale
    
@patch_to(CKKSEncoder)
def encode(self, z: np.array) -> Polynomial:
    """Encodes a vector by expanding it first to H,
    scale it, project it on the lattice of sigma(R), and performs
    sigma inverse.
    """
    pi_z = self.pi_inverse(z)
    print("pi_z:",pi_z)
    scaled_pi_z = self.scale * pi_z
    print("scaled_pi_z:",scaled_pi_z)
    rounded_scale_pi_zi = self.sigma_R_discretization(scaled_pi_z)
    print("rounded_scale_pi_zi:",rounded_scale_pi_zi)
    p = self.sigma_inverse(rounded_scale_pi_zi)
    print("p:",p)
    # We round it afterwards due to numerical imprecision
    coef = np.round(np.real(p.coef)).astype(int)
    print("coef:",coef)
    p = Polynomial(coef)
    print("p:",p)
    return p

@patch_to(CKKSEncoder)
def decode(self, p: Polynomial) -> np.array:
    """Decodes a polynomial by removing the scale, 
    evaluating on the roots, and project it on C^(N/2)"""
    rescaled_p = p / self.scale
    print("rescaled_p:",rescaled_p)
    z = self.sigma(rescaled_p)
    print("z:",z)
    pi_z = self.pi(z)
    return pi_z

scale = 64
M=8
encoder = CKKSEncoder(M, scale)

#测试
z = np.array([3 +4j, 2 - 1j])
print("明文：",z)
p = encoder.encode(z)
print("编码后：",p)
d = encoder.decode(p)
print("解码后：",d)
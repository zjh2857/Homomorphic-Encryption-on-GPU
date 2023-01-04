primes = []
prime = 2 ** 40

for i in range(8):
    prime = next_prime(prime)
    primes.append(prime)
print(primes)
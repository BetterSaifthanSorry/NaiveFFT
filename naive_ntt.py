# ntt operates much like fft 
# it operates over an integer "feild" which just mean integers modulo some large prime p
# and the multiplied polynomials live inside a ring

n = 32

def findInv(n):
    for i in range(0, q):
        if (n * i) % q == 1:
            return i

def modSquare(k, q):
    return k**n % q

def prime(n):
    for i in range(2, round(n ** 0.5)):
        if n % i == 0:
            return False
    return True

def nearestPowerOf2(n):
    prod = 1
    while(prod < n):
        prod = prod * 2
    return prod

# for the math to work out, q has to be such that n | (q - 1) and prime
# we can pick a q such that x*n + 1
# lets find it through code

def getq(n):
    for x in range(n + 1, 2**32):
        q = x * n + 1
        if (prime(q) and (q-1) % n == 0):
            return q
    return -1 # doesn't exist

q = getq(n)

ntt_roots = [
    1, 722, 282, 761, 753, 660, 449, 116, 
    468, 1131, 739, 752, 720, 903, 569, 252, 
    1152, 431, 871, 392, 400, 493, 704, 1037, 
    685, 22, 414, 401, 433, 250, 584, 901
]

ntt_roots = [1, 286, 1086, 439, 1030, 565, 170, 194, 140, 838, 997, 351, 75, 696, 740, 641, 1152, 867, 67, 714, 123, 588, 983, 959, 1013, 315, 156, 802, 1078, 457, 413, 512]
print(f"n: {n} q: {q}")

def ntt(polynomialVec):
    if len(polynomialVec) == 1:
        return [polynomialVec[0]]
    n = len(polynomialVec)
    
    evenCoeffs = [polynomialVec[i] for i in range(0, n, 2)]
    oddCoeffs = [polynomialVec[i] for i in range(1, n, 2)]

    odd, even = ntt(oddCoeffs), ntt(evenCoeffs)

    k = n // 2
    y = [0] * n
    roots = [ntt_roots[i] for i in range(0, len(ntt_roots), len(ntt_roots) // n)]
    for i in range(0, k): 
        y[i] = even[i] + roots[i]*odd[i] % q
        y[i+k] = even[i] + roots[i+k]*odd[i] % q
    return y

def intt(sampleVec):
    if len(sampleVec) == 1:
        return [sampleVec[0]]
    n = len(sampleVec)
    
    evenCoeffs = [sampleVec[i] for i in range(0, n, 2)]
    oddCoeffs = [sampleVec[i] for i in range(1, n, 2)]

    odd, even = intt(oddCoeffs), intt(evenCoeffs)

    k = n // 2
    y = [0] * n
    roots = [pow(ntt_roots[i], -1, q) for i in range(0, len(ntt_roots), len(ntt_roots) // n)]
    for i in range(0, k): 
        y[i] = even[i] + roots[i]*odd[i] % q
        y[i+k] = even[i] + roots[i+k]*odd[i] % q
    return y

def pad(vec, n):
    N = len(vec)
    return vec + ([0] * (n-N))

def fastPolynomialMultiply(polynomialA, polynomialB):
    nA, nB = len(polynomialA), len(polynomialB)
    n = nearestPowerOf2(nA + nB - 1)
    polynomialA, polynomialB = pad(polynomialA, n), pad(polynomialB, n)
    fft_a, fft_b = ntt(polynomialA), ntt(polynomialB)
    n = len(fft_a)
    fft_ab = []
    for (sample_a, sample_b) in zip(fft_a, fft_b):
        fft_ab.append(sample_a * sample_b)
    polynomialAB = intt(fft_ab)
    n_inv = findInv(n)
    polynomialAB = [(x * n_inv) % q for x in polynomialAB]
    return polynomialAB

def posWrappedPolynomialMultiply(polynomialA, polynomialB):
    nA, nB = len(polynomialA), len(polynomialB)
    n = nearestPowerOf2(max(nA, nB) + 1)
    polynomialA, polynomialB = pad(polynomialA, n), pad(polynomialB, n)
    fft_a, fft_b = ntt(polynomialA), ntt(polynomialB)
    n = len(fft_a)
    fft_ab = []
    for (sample_a, sample_b) in zip(fft_a, fft_b):
        fft_ab.append(sample_a * sample_b)
    polynomialAB = intt(fft_ab)
    n_inv = findInv(n)
    polynomialAB = [(x * n_inv) % q for x in polynomialAB]
    return polynomialAB

psi_table = [1, 943, 286, 1049, 1086, 234, 439, 50, 1030, 464, 565, 109, 170, 43, 194, 768, 140, 578, 838, 429, 997, 476, 351, 82, 75, 392, 696, 271, 740, 255, 641, 291]
psi_inv_table = [1, 862, 512, 898, 413, 882, 457, 761, 1078, 1071, 802, 677, 156, 724, 315, 575, 1013, 385, 959, 1110, 983, 1044, 588, 689, 123, 1103, 714, 919, 67, 104, 867, 210]

def negacyclicPolynomialMultiply(polynomialA, polynomialB):
    for i in range(len(polynomialA)):
        polynomialA[i] = (polynomialA[i] * psi_table[i*8]) % q
    for i in range(len(polynomialB)):
        polynomialB[i] = (polynomialB[i] * psi_table[i*8]) % q
    polynomialAB = posWrappedPolynomialMultiply(polynomialA, polynomialB)

    for i in range(len(polynomialAB)):
        polynomialAB[i] = (polynomialAB[i] * psi_inv_table[i*8]) % q
    return polynomialAB
    

print(fastPolynomialMultiply([1, 2, 3], [4, 6, 8]))
print(f"Positive: {posWrappedPolynomialMultiply([1, 2, 3], [4, 6, 8])}")
print(f"Negacyclic: {negacyclicPolynomialMultiply([1, 2, 3], [4, 6, 8])}")

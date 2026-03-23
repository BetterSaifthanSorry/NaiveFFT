import random
import cmath

degree = int(input("Enter the desired degree of the polynomial (Coefficent Vector View): "))

polynomialA = [random.randint(0, 100) for _ in range(degree + 1)]
polynomialB = [random.randint(0, 100) for _ in range(degree + 1)]

def printPolynomial(polynomialVector):
    variables = [""]
    for exp in range(1, len(polynomialVector) + 1):
        variables.append(f"x^{exp}")
    polynomial = [f"{coeff}{var}" for (coeff, var) in zip(polynomialVector, variables)]
    return " + ".join(polynomial)

print(f"Polynomial A: {printPolynomial(polynomialA)}")
print(f"Polynomial B: {printPolynomial(polynomialB)}")

# this is the naive way to multiply two polynomials in coefficent vector form
# this is quadratic time - O(n^2)
def naiveVectorMultiply(polynomialA, polynomialB):
    # schoolbook convolution
    # y[k] = (g*h)[k] += g[i]h[k-i] for i = 0 to i = k
    n = len(polynomialA)
    degree = n - 1
    # pad with zeroes
    polynomialA = polynomialA + ([0] * degree)
    polynomialB = polynomialB + ([0] * degree)
    # res = [0] * (2 * degree + 1)
    res = []

    for deg in range(0, 2*degree + 1):
        prod = 0
        for i in range(0, deg + 1):
            prod += polynomialA[i] * polynomialB[deg - i]
        res.append(prod)
    return res

print(f"Multiplied Polynomial: {printPolynomial(naiveVectorMultiply(polynomialA, polynomialB))}")
    
# let's take a little detour now and introduce another view 
# the sample view

# to take the sample view of a polynomial f(x) of degree d you need d+1 samples of the form (x, f(x))
# any polynomial of degree d is uniquely characterized by d samples

degree = int(input("Enter the desired degree of the polynomial (Sample View): "))

polynomialA = [random.randint(0, 100) for _ in range(degree + 1)]
polynomialB = [random.randint(0, 100) for _ in range(degree + 1)]

# let us now evaluate these polynomials at d + 1 random values

def evaluate(polynomialVector, x):
    degree = len(polynomialVector) -  1
    powersOfX = [x**i for i in range(0, degree + 1)]
    # now we can just take the dot product of these two vectors
    res = 0
    for (coeff, var) in zip(polynomialVector, powersOfX):
        res += coeff * var
    return res

samplePoints = [random.randint(0, 100) for _ in range(degree + 1)]

samplesA = [(i, evaluate(polynomialA, i)) for i in samplePoints]
samplesB = [(i, evaluate(polynomialB, i)) for i in samplePoints]

print(f"Sample View of Polynomial A: {samplesA}")
print(f"Sample View of Polynomial B: {samplesB}")

# now lets try to multiply two polynomials in the sample universe

def sampleMultiply(samplesA, samplesB):
    return [(x, A*B) for ((x, A),(_, B)) in zip(samplesA, samplesB)]

productSamples = sampleMultiply(samplesA, samplesB)
print(f"Sample View of Product: {productSamples}")

# let's confirm the product of samples is equal to the samples of the product

prodCoefficents = naiveVectorMultiply(polynomialA, polynomialB)
samplesProd = [(i, evaluate(prodCoefficents, i)) for i in samplePoints]
print(productSamples == samplesProd)

# so multiplication is quite cheaper in the sample world...
# it's O(n)

# so if we convert from the vector to the sample view we can get multiplication for cheaper...
# but if you scroll up, you'll see what that conversion entails
# evaluating a polynomial of degree n-1 at n points costs O(n^2) 

# so the O(n^2) conversion dominates the O(n) multiplication

# what if we could get the evaluations for cheaper...
# the insight comes from evaluating at special points
# points where you could reuse certain values...
# those special points are the roots of unity...

def getRootsOfUnity(n):
    roots = []
    for k in range(n):
        # 2j is Python's way of writing 2 * sqrt(-1)
        angle = 2 * cmath.pi * k / n
        root = cmath.exp(angle * 1j) 
        roots.append(root)
    return roots

def getConjugateRootsOfUnity(n):
    roots = []
    for k in range(n):
        # 2j is Python's way of writing 2 * sqrt(-1)
        angle = 2 * cmath.pi * k / n
        root = cmath.exp(-angle * 1j) 
        roots.append(root)
    return roots

# print(f"Roots of Unity: {getRootsOfUnity(32)}")

# the special nature of the roots of unity allows us to do the conversion in O(nlogn)
# through a special divide-and-conquer algorithm called the Fast-Fourier Transform

def nearestPowerOf2(n):
    prod = 1
    while(prod < n):
        prod = prod * 2
    return prod



def fft(polynomialVec):
    if len(polynomialVec) == 1:
        return [polynomialVec[0]]
    n = len(polynomialVec)
    
    evenCoeffs = [polynomialVec[i] for i in range(0, n, 2)]
    oddCoeffs = [polynomialVec[i] for i in range(1, n, 2)]

    odd, even = fft(oddCoeffs), fft(evenCoeffs)

    k = n // 2
    y = [0] * n
    roots = getRootsOfUnity(n)
    for i in range(0, k): 
        y[i] = even[i] + roots[i]*odd[i]
        y[i+k] = even[i] + roots[i+k]*odd[i]
    return y

def ifft(sampleVec):
    if len(sampleVec) == 1:
        return [sampleVec[0]]
    n = len(sampleVec)
    
    evenCoeffs = [sampleVec[i] for i in range(0, n, 2)]
    oddCoeffs = [sampleVec[i] for i in range(1, n, 2)]

    odd, even = ifft(oddCoeffs), ifft(evenCoeffs)

    k = n // 2
    y = [0] * n
    roots = getConjugateRootsOfUnity(n)
    for i in range(0, k): 
        y[i] = even[i] + roots[i]*odd[i]
        y[i+k] = even[i] + roots[i+k]*odd[i]
    return y

def pad(vec, n):
    N = len(vec)
    return vec + ([0] * (n-N))

# let us now write the final O(nlogn) polynomial multiplication routine using fft as a subroutine
def fastPolynomialMultiply(polynomialA, polynomialB):
    nA, nB = len(polynomialA), len(polynomialB)
    finalDegree = (nA - 1) + (nB - 1)
    n = nearestPowerOf2(finalDegree + 1)
    polynomialA, polynomialB = pad(polynomialA, n), pad(polynomialB, n)
    print(f"Polynomial A: {polynomialA}")
    print(f"Polynomial B: {polynomialB}")
    fft_a, fft_b = fft(polynomialA), fft(polynomialB)
    n = len(fft_a)
    fft_ab = []
    for (sample_a, sample_b) in zip(fft_a, fft_b):
        fft_ab.append(sample_a * sample_b)
    polynomialAB = ifft(fft_ab)
    polynomialAB = [round((x/n).real) for x in polynomialAB]
    return polynomialAB[:finalDegree+1]

print(fastPolynomialMultiply([1, 2, 3], [4, 6, 8]))
print(naiveVectorMultiply([1, 2, 3], [4, 6, 8]))

print(naiveVectorMultiply([1, 2, 3], [4, 6, 8]) == fastPolynomialMultiply([1, 2, 3], [4, 6, 8]))
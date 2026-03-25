#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

const int n = 32;
int q;

int ntt_roots[32] = {1, 286, 1086, 439, 1030, 565, 170, 194, 140, 838, 997, 351, 75, 696, 740, 641, 1152, 867, 67, 714, 123, 588, 983, 959, 1013, 315, 156, 802, 1078, 457, 413, 512};

int find_inverse_mod_q(int n, int q)
{
    for (int i = 0; i < q; i++)
    {
        if ((n * i) % q == 1)
            return i;
    }
    return -1;
}

bool is_prime(int n)
{
    for (int i = 2; i < sqrt(n); i++)
    {
        if (n % i == 0)
            return false;
    }
    return true;
}

int find_q(int n)
{
    for (unsigned int k = n + 1; k < (unsigned int)-1; k++)
    {
        q = k * n + 1;
        if (is_prime(q) && ((q - 1) % n == 0))
            return q;
    }
    return -1;
}

int find_inv(int n)
{
    for (int i = 0; i < q; i++)
    {
        if ((n * i) % q == 1)
            return i;
    }
    return -1;
}

int nearest_power_of_2(int n)
{
    int prod = 1;
    while (prod < n)
        prod *= 2;
    return prod;
}

int pow_int(int n, int b)
{
    int prod = 1;
    for (int i = 0; i < b; i++)
        prod *= n;
    return prod;
}

int mod_pow_int(int n, int b, int q)
{
    int prod = 1;
    for (int i = 0; i < b; i++)
        prod *= n;
    return prod % q;
}

int to_binary(int n, int *arr)
{
    while (n != 0)
    {
        printf("n: %d ", n);
        int i = 0;
        while (pow_int(2, i) <= n)
            i++;
        --i;
        printf("setting i to %d\n", i);
        arr[i] = 1;
        n = n % pow_int(2, i);
    }
}

int from_binary(int *arr, int len)
{
    int res = 0;
    for (int i = 0; i < len; i++)
    {
        if (arr[i] != 0)
            res += pow_int(2, i);
    }
    return res;
}

int invert_binary(int n)
{
    int arr[32];
    to_binary(n, arr);
    for (int i = 0; i < 16; i++)
    {
        int temp = arr[i];
        arr[i] = arr[i + 16];
        arr[i + 16] = temp;
    }
    int rev = from_binary(arr, 32);
    return rev;
}

// in C the in place Cooley Tukey is easier

int *NTT(int *polynomialVec, int len, int padded_len)
{
    int *padded_polynomial = (int *)malloc(padded_len * sizeof(int));
    memset(padded_polynomial, 0, padded_len * sizeof(int));
    for (int i = 0; i < len; i++) padded_polynomial[i] = polynomialVec[i];
    for (int i = 1, j = 0; i < padded_len; i++)
    {
        int bit = padded_len >> 1;
        for (; j & bit; bit >>= 1)
        {
            j ^= bit;
        }
        j ^= bit;

        if (i < j)
        {
            int temp = padded_polynomial[i];
            padded_polynomial[i] = padded_polynomial[j];
            padded_polynomial[j] = temp;
        }
    }

    // cooley-tukey in place

    int stride = 2;
    while (stride <= padded_len)
    {
        for (int i = 0; i < padded_len; i += stride)
        {
            int l = 32 / stride;
            for (int j = 0; j < stride / 2; j++)
            {
                int a = i + j;
                int b = i + j + stride / 2;

                long long original_A = padded_polynomial[a];
                long long mixed_term = (1LL * ntt_roots[j * l] * padded_polynomial[b]) % q;

                padded_polynomial[a] = (original_A + mixed_term) % q;
                padded_polynomial[b] = (original_A - mixed_term + q) % q;
            }
        }
        stride *= 2;
    }
    return padded_polynomial;
}

int *INTT(int *polynomialVec, int len, int padded_len)
{
    int *padded_polynomial = (int *)malloc(padded_len * sizeof(int));
    memset(padded_polynomial, 0, padded_len * sizeof(int));
    for (int i = 0; i < len; i++) padded_polynomial[i] = polynomialVec[i];
    for (int i = 1, j = 0; i < padded_len; i++)
    {
        int bit = padded_len >> 1;
        for (; j & bit; bit >>= 1)
        {
            j ^= bit;
        }
        j ^= bit;

        if (i < j)
        {
            int temp = padded_polynomial[i];
            padded_polynomial[i] = padded_polynomial[j];
            padded_polynomial[j] = temp;
        }
    }

    // cooley-tukey in place

    int stride = 2;
    while (stride <= padded_len)
    {
        for (int i = 0; i < padded_len; i += stride)
        {
            int l = 32 / stride;
            for (int j = 0; j < stride / 2; j++)
            {
                int a = i + j;
                int b = i + j + stride / 2;

                long long original_A = padded_polynomial[a];
                long long mixed_term = (1LL * find_inv(ntt_roots[j * l]) * padded_polynomial[b]) % q;

                padded_polynomial[a] = (original_A + mixed_term) % q;
                padded_polynomial[b] = (original_A - mixed_term + q) % q;
            }
        }
        stride *= 2;
    }
    int inv_n = find_inv(len);
    for (int i = 0; i< len; i++)
        padded_polynomial[i] = (padded_polynomial[i] * inv_n) % q; 
    return padded_polynomial;
}

void print_array(int *arr, int len)
{
    for (int i = 0; i < len; i++)
        printf("%d, ", arr[i]);
    printf("\n");
}

int *positivelyWrappedNTT(int *polynomial_a, int *polynomial_b)
{
    
    int *ntt_a = NTT(polynomial_a, 4, 8);
    int *ntt_b = NTT(polynomial_b, 4, 8);

    for (int i = 0; i < 8; i++)
    {
        ntt_a[i] =(1LL * ntt_a[i] * ntt_b[i]) % q;
    }

    int *c = INTT(ntt_a, 8, 8);

    return c;
}

int *fastPolynomialMultiply(int *polynomial_a, int len_a, int *polynomial_b, int len_b)
{
    int deg_a = len_a - 1;
    int deg_b = len_b - 1;
    int deg_c = deg_a + deg_b;
    int len_c = nearest_power_of_2(deg_c + 1);
    int *ntt_a = NTT(polynomial_a, len_a, nearest_power_of_2(len_c));
    int *ntt_b = NTT(polynomial_b, len_b, nearest_power_of_2(len_c));

    for (int i = 0; i < len_c; i++)
    {
        ntt_a[i] =(1LL * ntt_a[i] * ntt_b[i]) % q;
    }

    int *c = INTT(ntt_a, len_c, len_c);

    return c;
}

int main()
{
    q = find_q(n);
    int arr[32] = {0};
    int n = 32;
    to_binary(n, arr);
    int a[4] = {1, 2, 3, 0};
    int b[4] = {4, 6, 8, 0};
    print_array(positivelyWrappedNTT(a, b), 4);
    print_array(fastPolynomialMultiply(a, 4, b, 4), 8);
}
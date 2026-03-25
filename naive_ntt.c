#include <stdio.h>
#include <math.h>
#include <stdbool.h>

const int n = 32;
int q;

int ntt_roots[n] = {1, 286, 1086, 439, 1030, 565, 170, 194, 140, 838, 997, 351, 75, 696, 740, 641, 1152, 867, 67, 714, 123, 588, 983, 959, 1013, 315, 156, 802, 1078, 457, 413, 512};

int find_inverse_mod_q(int n, int q){
    for (int i = 0; i < q; i++){
        if ((n * i) % q == 1)
            return i;
    }
    return -1;
}

bool is_prime(int n){
    for (int i = 2; i < sqrt(n); i++){
        if (n % i == 0)
            return false;
    }
    return true;
}

int find_q(int n){
    for (unsigned int k = n + 1; k < (unsigned int)-1; k++){
        q = k * n + 1;
        if (is_prime(q) && ((q-1) % n == 0))
            return q;
    }
    return -1;
}

int find_inv(int n){
    for (int i = 0; i < q; i++){
        if ((n * i) % q == 1)
            return i;
    }
    return -1;
}

int nearest_power_of_2(int n){
    int prod = 1;
    while (prod < n)
        prod *= 2;
    return prod;
}

int NTT(int* polynomialVec, int len){
    int new_size = nearest_power_of_2(len);
    int padded_polynomial[new_size];
    
}

int main(){
    q = find_q(n);
    printf("n: %d q: %d", n, q);
}
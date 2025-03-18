#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

// Structure to represent an elliptic curve point (x, y)
typedef struct {
    mpz_t x;
    mpz_t y;
} point_t;

// Initialize a point
void init_point(point_t *p) {
    mpz_init(p->x);
    mpz_init(p->y);
}

// Clear a point
void clear_point(point_t *p) {
    mpz_clear(p->x);
    mpz_clear(p->y);
}

// Point addition: r = p + q on curve y^2 = x^3 + ax + b mod n
void point_add(point_t *r, point_t *p, point_t *q, mpz_t a, mpz_t n, int *factor_found, mpz_t factor) {
    mpz_t lambda, tmp, den;
    mpz_init(lambda);
    mpz_init(tmp);
    mpz_init(den);

    // If p == q, use point doubling
    if (mpz_cmp(p->x, q->x) == 0 && mpz_cmp(p->y, q->y) == 0) {
        mpz_set(den, p->y);
        mpz_mul_ui(den, den, 2); // den = 2y
        mpz_mul(tmp, p->x, p->x);
        mpz_mul_ui(tmp, tmp, 3); // tmp = 3x^2
        mpz_add(tmp, tmp, a);    // tmp = 3x^2 + a
    } else {
        mpz_sub(den, q->x, p->x); // den = x_q - x_p
        mpz_sub(tmp, q->y, p->y); // tmp = y_q - y_p
    }

    // Compute lambda = tmp / den mod n
    if (mpz_invert(lambda, den, n) == 0) {
        // Inverse doesn’t exist: factor found
        mpz_gcd(factor, den, n);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, n) < 0) {
            *factor_found = 1;
            printf("Found factor: ");
            mpz_out_str(stdout, 10, factor);
            printf("\n");
            fflush(stdout); // Ensure immediate output
        }
        goto cleanup;
    }
    mpz_mul(lambda, tmp, lambda);
    mpz_mod(lambda, lambda, n);

    // r_x = lambda^2 - p_x - q_x mod n
    mpz_mul(tmp, lambda, lambda);
    mpz_sub(tmp, tmp, p->x);
    mpz_sub(tmp, tmp, q->x);
    mpz_mod(r->x, tmp, n);

    // r_y = lambda * (p_x - r_x) - p_y mod n
    mpz_sub(tmp, p->x, r->x);
    mpz_mul(tmp, lambda, tmp);
    mpz_sub(tmp, tmp, p->y);
    mpz_mod(r->y, tmp, n);

cleanup:
    mpz_clear(lambda);
    mpz_clear(tmp);
    mpz_clear(den);
}

// Scalar multiplication: r = k * p using double-and-add method
void point_multiply(point_t *r, point_t *p, mpz_t k, mpz_t a, mpz_t n, int *factor_found, mpz_t factor) {
    point_t temp, result;
    init_point(&temp);
    init_point(&result);
    mpz_set(result.x, p->x);
    mpz_set(result.y, p->y);

    size_t bits = mpz_sizeinbase(k, 2);
    for (size_t i = 0; i < bits && !*factor_found; i++) {
        if (mpz_tstbit(k, i)) {
            point_add(&temp, &result, p, a, n, factor_found, factor);
            if (!*factor_found) {
                mpz_set(result.x, temp.x);
                mpz_set(result.y, temp.y);
            }
        }
        if (!*factor_found) {
            point_add(&temp, &result, &result, a, n, factor_found, factor);
            if (!*factor_found) {
                mpz_set(result.x, temp.x);
                mpz_set(result.y, temp.y);
            }
        }
    }
    if (!*factor_found) {
        mpz_set(r->x, result.x);
        mpz_set(r->y, result.y);
    }

    clear_point(&temp);
    clear_point(&result);
}

// ECM implementation with continuous factor checking
void ecm_factor(mpz_t n, unsigned long B1) {
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    point_t P, Q;
    mpz_t a, k, tmp, factor;
    init_point(&P);
    init_point(&Q);
    mpz_init(a);
    mpz_init(k);
    mpz_init(tmp);
    mpz_init(factor);

    int factor_found = 0;
    int attempts = 0;
    const int max_attempts = 100; // Increased to keep searching

    while (attempts < max_attempts) {
        // Choose random point and curve parameter a
        mpz_urandomm(P.x, state, n);
        mpz_urandomm(P.y, state, n);
        mpz_urandomm(a, state, n);

        // Compute k = product of small primes up to B1
        mpz_set_ui(k, 1);
        for (unsigned long p = 2; p <= B1; p++) {
            if (is_prime(p)) {
                unsigned long e = 1;
                unsigned long power = p;
                while (power <= B1) {
                    power *= p;
                    e++;
                }
                power /= p;
                mpz_mul_ui(tmp, k, power);
                mpz_set(k, tmp);
            }
        }

        // Perform scalar multiplication: Q = k * P
        factor_found = 0;
        point_multiply(&Q, &P, k, a, n, &factor_found, factor);

        if (factor_found) {
            // Divide n by the found factor and continue
            mpz_t quotient;
            mpz_init(quotient);
            mpz_divexact(quotient, n, factor);
            mpz_set(n, quotient);
            mpz_clear(quotient);
            printf("Remaining number: ");
            mpz_out_str(stdout, 10, n);
            printf("\n");
            fflush(stdout);
            if (mpz_cmp_ui(n, 1) == 0) {
                printf("Number fully factored.\n");
                break;
            }
        }

        attempts++;
        printf("Attempt %d completed.\n", attempts);
        fflush(stdout);
    }

    if (attempts >= max_attempts) {
        printf("No more factors found after %d attempts.\n", max_attempts);
    }

    clear_point(&P);
    clear_point(&Q);
    mpz_clear(a);
    mpz_clear(k);
    mpz_clear(tmp);
    mpz_clear(factor);
    gmp_randclear(state);
}

// Simple primality test for small numbers
int is_prime(unsigned long n) {
    if (n < 2) return 0;
    if (n == 2) return 1;
    if (n % 2 == 0) return 0;
    for (unsigned long i = 3; i * i <= n; i += 2) {
        if (n % i == 0) return 0;
    }
    return 1;
}

int main() {
    mpz_t n;
    mpz_init_set_str(n, "17135288265189914517928350174415883458808762336089003353130797688780057659774311825837762434272848522107627891868701632769900177103123828589399992283229465278080016870350078985858411427955061552777518854837895495586635215995279670451965754378437758118675571133179996475365825296232680234341551427580200803296460434550371268034913769321682298406519398630780956316496390694222865119928015212532501", 10);

    unsigned long B1 = 10000; // Larger B1 for bigger factors

    printf("Attempting to factor: ");
    mpz_out_str(stdout, 10, n);
    printf("\n");
    fflush(stdout);

    ecm_factor(n, B1);

    mpz_clear(n);
    return 0;
}
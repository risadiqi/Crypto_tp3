//
//  TP6_RSA
//  

#include <stdio.h>
#include <iostream>
#include <gmp.h>
#include <time.h>

#define BITSTRENGTH  14              /* size of modulus (n) in bits */
#define PRIMESIZE (BITSTRENGTH / 2)  /* size of the primes p and q  */

/* Declare global variables */

mpz_t d,e,n;
mpz_t M,C;


void Expo_By_Squaring(mpz_t result, mpz_t g, mpz_t k, mpz_t p){

    if(mpz_sgn(k) < 0){
        mpz_invert(g, g, p);
        mpz_neg(k, k);
    }

    else if(mpz_sgn(k) == 0){
        mpz_set_ui(result, (unsigned long int)1);    
    }

    else {
        mpz_t y;
        mpz_init_set_ui(y, (unsigned long int)1);

        while(mpz_cmp_ui(k, 1) > 0){

            if(mpz_even_p(k) != 0){
                mpz_powm_ui(g, g, (unsigned long int)2, p);
                mpz_divexact_ui(k, k, 2);
            }

            else {
                mpz_mul(y, g, y);
                mpz_mul(g, g, g);
                mpz_sub_ui(k, k, (unsigned long int)1);
                mpz_divexact_ui(k, k, 2);
            }
        }
        mpz_mul(result, g, y);
        mpz_mod(result, result, p);
    }

}

/* Main subroutine */
int main()
{
    /* Initialize the GMP integers */
    mpz_init(d);
    mpz_init(e);
    mpz_init(n);
    
    mpz_init(M);
    mpz_init(C);
    
 
    /* This function creates the keys. The basic algorithm is...
     *
     *  1. Generate two large distinct primes p and q randomly
     *  2. Calculate n = pq and x = (p-1)(q-1)
     *  3. Select a random integer e (1<e<x) such that gcd(e,x) = 1
     *  4. Calculate the unique d such that ed = 1(mod x)
     *  5. Public key pair : (e,n), Private key pair : (d,n)
     *
     */
    
    /*
     *  Step 1 : Get two large primes.
     */
    mpz_t p,q;
    mpz_t randp, randq;
    mpz_init(randp);
    mpz_init(randq);
    mpz_init(p);
    mpz_init(q);

    gmp_randstate_t r;
    gmp_randinit_default(r);
    gmp_randseed_ui(r, time(NULL));

    mpz_urandomb(randp, r, PRIMESIZE);
    mpz_urandomb(randq, r, PRIMESIZE);

    mpz_nextprime(p, randp);
    mpz_nextprime(q, randq);
    
    char p_str[1000];
    char q_str[1000];
    mpz_get_str(p_str,10,p);
    mpz_get_str(q_str,10,q);
    
    std::cout << "Random Prime 'p' = " << p_str <<  std::endl;
    std::cout << "Random Prime 'q' = " << q_str <<  std::endl;
    
    /*
     *  Step 2 : Calculate n (=pq) ie the 1024 bit modulus
     *  and x (=(p-1)(q-1)).
     */
    char n_str[1000];
    mpz_t x;
    mpz_init(x);

    /* Calculate n... */
    mpz_mul(n,p,q);
    mpz_get_str(n_str,10,n);
    std::cout << "\t n = " << n_str << std::endl;
    
    
    /* Calculate x... */
    mpz_t p_minus_1,q_minus_1;
    mpz_init(p_minus_1);
    mpz_init(q_minus_1);

    mpz_sub_ui(p_minus_1,p,(unsigned long int)1);
    mpz_sub_ui(q_minus_1,q,(unsigned long int)1);

    mpz_mul(x,p_minus_1,q_minus_1);
    char phi_str[1000];
    mpz_get_str(phi_str,10,x);
    std::cout << "\t phi(n) = " << phi_str << std::endl;

    /*
     *  Step 3 : Get small odd integer e such that gcd(e,x) = 1.
     */
    mpz_init(e);
    mpz_t rop;

    mpz_init(rop);

    do{
        gmp_randstate_t ra;
        gmp_randinit_default(ra);
        gmp_randseed_ui(ra, time(NULL));

        mpz_urandomb(e, ra, PRIMESIZE);
        mpz_gcd(rop, e, x);

    } while (mpz_cmp_ui(rop, 1) != 0);

    char e_str[1000];
    mpz_get_str(e_str,10,e);
    std::cout << "\t e = " << e_str << std::endl;

    /*
     *  Step 4 : Calculate unique d such that ed = 1(mod x)
     */
    mpz_init(d);
    mpz_invert(d, e, x);
    
    char d_str[1000];
    mpz_get_str(d_str,10,d);
    std::cout << "\t d = " << d_str << std::endl << std::endl;

    /*
     *  Step 5 : Print the public and private key pairs...
     */
    std::cout << "Public Keys  (e,n): ( " << e_str <<" , " << n_str << " )" << std::endl;
    std::cout << "Private Keys (d,n): ( " << d_str <<" , " << n_str << " )" << std::endl;
    /*
     *  Encrypt
     */

    /* 
     * Step 6: Generate random message M < n 
     */
    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, time(NULL));

    mpz_urandomb(M, rs, BITSTRENGTH - 1);  // M aléatoire de taille inférieure à n

    // Vérifier que M < n, sinon générer un nouveau
    while (mpz_cmp(M, n) >= 0) {
        mpz_urandomb(M, rs, BITSTRENGTH - 1);
    }

    char M_str[1000];
    mpz_get_str(M_str, 10, M);
    std::cout << "Message aléatoire M à chiffrer : " << M_str << std::endl;

    /*
     * Step 7: Encrypt the message M to get ciphertext C 
     */

    Expo_By_Squaring(C, M, e, n);  // C = M^e mod n

    char C_str[1000];
    mpz_get_str(C_str, 10, C);
    std::cout << "Message chiffré C : " << C_str << std::endl;

    /* 
     * Step 8: Decrypt the ciphertext C to get message M'
     */
    mpz_t M_decrypted;
    mpz_init(M_decrypted);
    Expo_By_Squaring(M_decrypted, C, d, n);  // M' = C^d mod n

    char M_decrypted_str[1000];
    mpz_get_str(M_decrypted_str, 10, M_decrypted);
    std::cout << "Message déchiffré M' : " << M_decrypted_str << std::endl;
    
    /* 
     * Clean up the GMP integers 
     */
    
    mpz_clear(p_minus_1);
    mpz_clear(q_minus_1);
    mpz_clear(x);
    mpz_clear(p);
    mpz_clear(q);
    
    mpz_clear(d);
    mpz_clear(e);
    mpz_clear(n);
    
    mpz_clear(M);
    mpz_clear(C);

    mpz_clear(rop);
    mpz_clear(randp);
    mpz_clear(randq);
}


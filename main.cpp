#include <stdio.h>   
#include <iostream>  
#include <gmp.h>     
#include <time.h>    

#define BITSTRENGTH  14              /* Taille de n en bits (modulus) */
#define PRIMESIZE (BITSTRENGTH / 2)  /* Taille des nombres premiers p et q en bits */

/* Déclaration des variables globales */
mpz_t d, e, n;  // d = exposant privé, e = exposant public, n = module
mpz_t M, C;     // M = message, C = message chiffré


// Fonction Expo_By_Squaring : Effectue l'exponentiation rapide pour calculer (g^k) mod p de manière efficace.
 
void Expo_By_Squaring(mpz_t result, mpz_t g, mpz_t k, mpz_t p) {
    if(mpz_sgn(k) < 0){  // Si l'exposant k est négatif, on prend l'inverse de g mod p
        mpz_invert(g, g, p);
        mpz_neg(k, k);
    }
    else if(mpz_sgn(k) == 0){  // Si k est nul, le résultat est 1
        mpz_set_ui(result, (unsigned long int)1);    
    }
    else {
        mpz_t y;
        mpz_init_set_ui(y, (unsigned long int)1);  // Initialisation de y = 1

        // Boucle pour l'exponentiation rapide
        while(mpz_cmp_ui(k, 1) > 0){
            if(mpz_even_p(k) != 0){  // Si k est pair
                mpz_mul(g, g, g);    // g = g^2
                mpz_mod(g, g, p);    // g = g mod p
                mpz_divexact_ui(k, k, 2);  // k = k / 2
            }
            else {  // Si k est impair
                mpz_mul(y, g, y);    // y = g * y
                mpz_mul(g, g, g);    // g = g^2
                mpz_sub_ui(k, k, (unsigned long int)1);  // k = k - 1
                mpz_divexact_ui(k, k, 2);  // k = k / 2
            }
        }
        mpz_mul(result, g, y);  // result = g * y
        mpz_mod(result, result, p);  // result = result mod p
    }
}


 //Fonction rabin_miller : implémente le test de primalité probabiliste de Rabin-Miller pour vérifier si n est premier.
 //'k' est le nombre d'itérations du test pour garantir une bonne précision.
 
bool rabin_miller(mpz_t n, int k) {
    if (mpz_even_p(n))  // Si n est pair, il n'est pas premier
        return false;

    mpz_t n_minus_1, t;
    mpz_init(n_minus_1);
    mpz_sub_ui(n_minus_1, n, 1);  // n_minus_1 = n - 1
    mpz_init_set(t, n_minus_1);  // Initialisation de t = n - 1

    int s = 0;
    // Divise n-1 par 2 jusqu'à ce qu'on obtienne un nombre impair
    while (mpz_even_p(t)) { 
        mpz_fdiv_q_ui(t, t, 2);  // t = t / 2
        s++;
    }

    gmp_randstate_t rstate;  // Initialise l'état pour les nombres aléatoires
    gmp_randinit_default(rstate);

    mpz_t a, x;
    mpz_init(a);
    mpz_init(x);

    for (int i = 0; i < k; i++) {
        do {
            gmp_randseed_ui(rstate, time(NULL));  // Initialise la graine pour les nombres aléatoires
            mpz_urandomm(a, rstate, n_minus_1);   // Génère a aléatoire dans [2, n-1]
        } while (mpz_cmp_ui(a, 2) < 0 || mpz_cmp(a, n_minus_1) >= 0);

        // Effectue l'exponentiation rapide x = a^t mod n
        Expo_By_Squaring(x, a, t, n);

        if (mpz_cmp_ui(x, 1) != 0 && mpz_cmp(x, n_minus_1) != 0) {
            bool composite = true;
            for (int r = 1; r < s; r++) {
                mpz_mul(x, x, x);    // x = x^2
                mpz_mod(x, x, n);    // x = x mod n

                if (mpz_cmp_ui(x, 1) == 0) {  // Si x est 1, n est composé
                    mpz_clear(a);
                    mpz_clear(x);
                    mpz_clear(n_minus_1);
                    mpz_clear(t);
                    gmp_randclear(rstate);
                    return false;
                }

                if (mpz_cmp(x, n_minus_1) == 0) {  // Si x est n-1, continue à tester
                    composite = false;
                    break;
                }
            }

            if (composite) {  // Si aucune condition n'est remplie, n est composé
                mpz_clear(a);
                mpz_clear(x);
                mpz_clear(n_minus_1);
                mpz_clear(t);
                gmp_randclear(rstate);
                return false;
            }
        }
    }

    // Si le test n'a pas échoué après k itérations, n est probablement premier
    mpz_clear(a);
    mpz_clear(x);
    mpz_clear(n_minus_1);
    mpz_clear(t);
    gmp_randclear(rstate);
    return true;
}

//Fonction gcd_euclidian : implémente l'algorithme d'Euclide pour calculer le PGCD de a et b.
 
void gcd_euclidian(mpz_t res, mpz_t a, mpz_t b) {
    if (mpz_cmp_ui(b, 0) == 0)
        mpz_set(res, a);  // Si b = 0, le PGCD est a
    else {
        mpz_t m;
        mpz_init(m);
        mpz_mod(m, a, b);  // m = a mod b
        gcd_euclidian(res, b, m);  // Appel récursif
    }
}


// Fonction principale : Génère les clés RSA, chiffre et déchiffre un message.

int main() {
    /* Initialisation des variables GMP */
    mpz_init(d);
    mpz_init(e);
    mpz_init(n);
    mpz_init(M);
    mpz_init(C);

    /* Étape 1 : Génération de deux grands nombres premiers p et q */
    mpz_t p, q;
    mpz_init(p);
    mpz_init(q);

    gmp_randstate_t r;
    gmp_randinit_default(r);
    gmp_randseed_ui(r, time(NULL));

    mpz_urandomb(p, r, PRIMESIZE);  // Génère un nombre premier p
    if (mpz_even_p(p)) {
        mpz_add_ui(p, p, 1);  // Assure que p est impair
    }

    mpz_urandomb(q, r, PRIMESIZE);  // Génère un nombre premier q
    if (mpz_even_p(q)) {
        mpz_add_ui(q, q, 1);  // Assure que q est impair
    }

    // Teste si p et q sont premiers avec Rabin-Miller
    while (rabin_miller(p, 7) == false) {
        mpz_add_ui(p, p, (unsigned long int)2);
    }

    while (rabin_miller(q, 7) == false || mpz_cmp(p, q) == 0) {
        mpz_add_ui(q, q, (unsigned long int)2);
    }

    /* Étape 2 : Calculer n = p * q et ϕ(n) = (p-1) * (q-1) */
    mpz_t x;
    mpz_init(x);
    mpz_mul(n, p, q);  // n = p * q

    mpz_t p_minus_1, q_minus_1;
    mpz_init(p_minus_1);
    mpz_init(q_minus_1);
    mpz_sub_ui(p_minus_1, p, 1);  // p-1
    mpz_sub_ui(q_minus_1, q, 1);  // q-1
    mpz_mul(x, p_minus_1, q_minus_1);  // x = (p-1) * (q-1)

    /* Étape 3 : Choisir e tel que gcd(e, x) = 1 */
    mpz_init(e);
    mpz_t rop;
    mpz_init(rop);
    do {
        gmp_randstate_t ra;
        gmp_randinit_default(ra);
        gmp_randseed_ui(ra, time(NULL));
        mpz_urandomb(e, ra, PRIMESIZE);  // Génère e
        gcd_euclidian(rop, e, x);  // Vérifie que gcd(e, x) = 1
    } while (mpz_cmp_ui(rop, 1) != 0);

    /* Étape 4 : Calculer d tel que e * d ≡ 1 (mod x) */
    mpz_invert(d, e, x);

    /* Étape 5 : Affichage des clés publiques et privées */
    std::cout << "Public Keys  (e,n): (" << mpz_get_str(NULL, 10, e) << ", " << mpz_get_str(NULL, 10, n) << ")" << std::endl;
    std::cout << "Private Keys (d,n): (" << mpz_get_str(NULL, 10, d) << ", " << mpz_get_str(NULL, 10, n) << ")" << std::endl;

    /* Étape 6 : Générer un message aléatoire M < n */
    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, time(NULL));
    mpz_urandomb(M, rs, BITSTRENGTH - 1);
    while (mpz_cmp(M, n) >= 0) {
        mpz_urandomb(M, rs, BITSTRENGTH - 1);  // S'assure que M < n
    }

    /* Étape 7 : Chiffrer le message M pour obtenir C */
    Expo_By_Squaring(C, M, e, n);  // C = M^e mod n

    /* Étape 8 : Déchiffrer le message chiffré C pour obtenir M' */
    mpz_t M_decrypted;
    mpz_init(M_decrypted);
    Expo_By_Squaring(M_decrypted, C, d, n);  // M' = C^d mod n

    /* Nettoyage des variables GMP */
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
    gmp_randclear(r);
}

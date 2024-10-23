#include <stdio.h>

#include "wreath/math/mpz.h"

int main(){
    wm_mpz_t v1;
    wm_mpz_init(v1);
    wm_mpz_ui(v1, 56);
    
    //wm_mpz_print_digits(v1);
    //printf("\n");
    wm_mpz_print(v1);
    printf("\n");
}

#include "wreath/math/mpz.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

void wm_mpz_realloc(wm_mpz_t z, size_t cap){
    static_assert(
        sizeof(wm_mpz_full_t)*CHAR_BIT/2 >= sizeof(wm_mpz_half_t)*CHAR_BIT,
        "Error: Invalid sizes for 'wm_mpz_full_t' and 'wm_mpz_half_t'"
    );
    static_assert(
        WM_MPZ_BASE_MASK > 10, 
        "Error: Invalid value for 'WM_MPZ_BASE_MASK' (must be greater than 10)"
    );

    if (cap < z->cap) return;
    wm_mpz_half_t* tmp = realloc(z->digits, sizeof(wm_mpz_half_t)*cap);
    if (!tmp) return;
    z->digits = tmp;
    z->cap = cap;
}
void wm_mpz_init(wm_mpz_t z){
    z->digits = NULL;
    z->sign = false;
    z->cap = 0;
    z->siz = 0;
    wm_mpz_realloc(z, WM_MPZ_MIN_DIGITS);
}
void wm_mpz_init_with_capacity(wm_mpz_t z, size_t cap){
    z->digits = NULL;
    z->sign = false;
    z->cap = 0;
    z->siz = 0;
    wm_mpz_realloc(z, cap);
}
void wm_mpz_clear(wm_mpz_t z){
    free(z->digits);
    z->digits = NULL;
}

void wm_mpz_set(wm_mpz_t out, wm_mpz_t rhs){
    wm_mpz_realloc(out, rhs->cap);
    memcpy(out->digits, rhs->digits, sizeof(wm_mpz_half_t)*out->cap);
    out->sign = rhs->sign;
    out->siz = rhs->siz;
}
void wm_mpz_ui(wm_mpz_t out, uintmax_t rhs){
    wm_mpz_realloc(out, WM_MPZ_MIN_DIGITS);
    out->sign = false;
    size_t a = 0;
    do{out->digits[a++] = rhs & WM_MPZ_BASE_MASK;}
    while (rhs >>= WM_MPZ_BASE_SHIFT);
    out->siz = a;
}
void wm_mpz_si(wm_mpz_t out, intmax_t rhs){
    wm_mpz_realloc(out, WM_MPZ_MIN_DIGITS);
    out->sign = rhs < 0;
    if (rhs < 0) rhs = -rhs; 
    size_t a = 0;
    do{out->digits[a++] = rhs & WM_MPZ_BASE_MASK;}
    while (rhs >>= WM_MPZ_BASE_SHIFT);
    out->siz = a;
}
void wm_mpz_swap(wm_mpz_t out, wm_mpz_t rhs){
    struct wm_mpz tmp;
    tmp.digits = out->digits;
    tmp.sign = out->sign;
    tmp.cap = out->cap;
    tmp.siz = out->siz;
    out->digits = rhs->digits;
    out->sign = rhs->sign;
    out->cap = rhs->cap;
    out->siz = rhs->siz;
    rhs->digits = tmp.digits;
    rhs->sign = tmp.sign;
    rhs->cap = tmp.cap;
    rhs->siz = rhs->siz;
}
int wm_mpz_cmp_raw(wm_mpz_t lhs, wm_mpz_t rhs){
    if (lhs->siz > rhs->siz) return 1;
    if (lhs->siz < rhs->siz) return -1;
    for (size_t a = lhs->siz-1; a < lhs->siz; a--){
        if (lhs->digits[a] > rhs->digits[a]) return 1;
        if (lhs->digits[a] < rhs->digits[a]) return -1;
    }
    return 0;
}

void wm_mpz_add_raw(wm_mpz_t out, wm_mpz_t lhs, wm_mpz_t rhs){
    if (lhs->siz < rhs->siz){
        wm_mpz_ptr tmp = lhs;
        lhs = rhs;
        rhs = tmp;
    }
    wm_mpz_realloc(out, lhs->siz+1);
    wm_mpz_half_t carry = 0;
    size_t a = 0;
    for (; a < rhs->siz; a++){
        carry += lhs->digits[a] + rhs->digits[a];
        out->digits[a] = carry & WM_MPZ_BASE_MASK;
        carry >>= WM_MPZ_BASE_SHIFT;
    }
    for (; a < lhs->siz; a++){
        carry += lhs->digits[a];
        out->digits[a] = carry & WM_MPZ_BASE_MASK;
        carry >>= WM_MPZ_BASE_SHIFT;
    }
    out->digits[a] = carry;
    out->siz = a+(carry!=0);
}
void wm_mpz_sub_raw(wm_mpz_t out, wm_mpz_t lhs, wm_mpz_t rhs){
    bool flip = false;
    if (wm_mpz_cmp_raw(lhs, rhs) == -1){
        wm_mpz_ptr tmp = lhs;
        lhs = rhs;
        rhs = tmp;
        flip = true;
    }
    wm_mpz_realloc(out, lhs->siz);
    wm_mpz_half_t borrow = 0;
    size_t a = 0;
    for (; a < rhs->siz; a++){
        borrow = lhs->digits[a] - rhs->digits[a] - borrow;
        out->digits[a] = borrow & WM_MPZ_BASE_MASK;
        borrow >>= WM_MPZ_BASE_SHIFT;
        borrow &= 1;
    }
    for (; borrow && a < lhs->siz; a++){
        borrow = lhs->digits[a] - borrow;
        out->digits[a] = borrow & WM_MPZ_BASE_MASK;
        borrow >>= WM_MPZ_BASE_SHIFT;
        borrow &= 1;
    }
    out->siz = a;
    for (; a < lhs->siz; a++) out->digits[a] = lhs->digits[a];
    out->sign = flip;
}

void wm_mpz_mul_raw_naive(wm_mpz_t out, wm_mpz_t lhs, wm_mpz_t rhs){
        wm_mpz_t out_tmp;
        wm_mpz_init_with_capacity(out_tmp, lhs->siz + rhs->siz);
        memset(out_tmp->digits, 0, sizeof(wm_mpz_half_t)*out_tmp->cap);
        wm_mpz_full_t carry = 0;
        wm_mpz_half_t* max_out = NULL;
        for (size_t a = 0; a < lhs->siz; a++){
            carry = 0;
            wm_mpz_full_t f = lhs->digits[a];
            wm_mpz_half_t* ptr_out = out_tmp->digits+a;
            wm_mpz_half_t* ptr_rhs = rhs->digits;
            wm_mpz_half_t* ptr_rhs_end = rhs->digits + rhs->siz;
            while (ptr_rhs < ptr_rhs_end){
                carry += *ptr_out + *(ptr_rhs++) * f;
                *(ptr_out++) = carry & WM_MPZ_BASE_MASK;
                carry >>= WM_MPZ_BASE_SHIFT;
            }
            if (carry) *ptr_out += carry & WM_MPZ_BASE_MASK;
            max_out = max_out > ptr_out ? max_out : ptr_out;
        }
        out_tmp->siz = (max_out - out_tmp->digits) + (carry != 0);
        wm_mpz_swap(out, out_tmp);
        wm_mpz_clear(out_tmp);
}
void wm_mpz_divrem_raw_naive(wm_mpz_t div, wm_mpz_t rem, wm_mpz_t lhs, wm_mpz_t rhs);
void wm_mpz_divrem_raw_single_naive(wm_mpz_t div, wm_mpz_half_t* rem, wm_mpz_t lhs, wm_mpz_half_t rhs){
    wm_mpz_t out_tmp;
    wm_mpz_init_with_capacity(out_tmp, lhs->siz);
    wm_mpz_full_t carry = 0;
    bool makes_smaller = rhs > lhs->digits[lhs->siz-1];
    if (makes_smaller && lhs->siz == 1){
        *rem = lhs->digits[0];
        wm_mpz_ui(div, 0);
        return;
    }
    if (makes_smaller) carry = lhs->digits[lhs->siz-1];
    out_tmp->siz = lhs->siz - makes_smaller;
    for (size_t a = lhs->siz-makes_smaller-1; a < lhs->siz; a--){
        wm_mpz_full_t tmp_div = lhs->digits[a] / rhs;
        carry <<= WM_MPZ_BASE_SHIFT;
        carry |= tmp_div;
        out_tmp->digits[a] = carry / rhs;
    }
    wm_mpz_swap(div, out_tmp);
    wm_mpz_clear(out_tmp);
    *rem = carry;
}

void wm_mpz_print_digits(wm_mpz_t z){
    if (z->sign) printf("- ");
    for (size_t a = 0; a < z->siz; a++){
        uintmax_t out = z->digits[a];
        printf("%jx ", out);
    }
}
void wm_mpz_print(wm_mpz_t z){
    wm_mpz_t div;
    wm_mpz_half_t rem;
    wm_mpz_init(div);
    wm_mpz_set(div, z);
    while (div->siz > 1 || div->digits[0] != 0){
        wm_mpz_divrem_raw_single_naive(div, &rem, div, 10);
        wm_mpz_print_digits(div);
        uintmax_t rem_tmp = rem;
        printf("%jx", rem_tmp);
    }
}

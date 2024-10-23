#ifndef WREATH_MATH_MPZ
#define WREATH_MATH_MPZ

#include <stdbool.h>
#include <stdint.h>

#ifndef wm_mpz_full_t
#define wm_mpz_full_t uint64_t
#endif
#ifndef wm_mpz_half_t
#define wm_mpz_half_t uint32_t
#endif

#ifndef WM_MPZ_BASE_SHIFT
#define WM_MPZ_BASE_SHIFT 10
#endif
#ifndef WM_MPZ_BASE_MASK
#define WM_MPZ_BASE_MASK (((wm_mpz_half_t)1 << WM_MPZ_BASE_SHIFT)-1)
#endif
#ifndef WM_MPZ_MIN_DIGITS
#define WM_MPZ_MIN_DIGITS 33
#endif

struct wm_mpz{
    wm_mpz_half_t* digits;
    size_t cap;
    size_t siz;
    bool sign;
};
typedef struct wm_mpz wm_mpz_t[1];
typedef struct wm_mpz* wm_mpz_ptr;
typedef const struct wm_mpz* wm_mpz_srcptr;

void wm_mpz_realloc(wm_mpz_t z, size_t cap);
void wm_mpz_init(wm_mpz_t z);
void wm_mpz_init_with_capacity(wm_mpz_t z, size_t cap);
void wm_mpz_clear(wm_mpz_t z);

void wm_mpz_set(wm_mpz_t out, wm_mpz_t rhs);
void wm_mpz_ui(wm_mpz_t out, uintmax_t rhs);
void wm_mpz_si(wm_mpz_t out, intmax_t rhs);
void wm_mpz_swap(wm_mpz_t out, wm_mpz_t rhs);

int wm_mpz_cmp_raw(wm_mpz_t lhs, wm_mpz_t rhs);

void wm_mpz_add_raw(wm_mpz_t out, wm_mpz_t lhs, wm_mpz_t rhs);
void wm_mpz_sub_raw(wm_mpz_t out, wm_mpz_t lhs, wm_mpz_t rhs);
void wm_mpz_mul_raw_naive(wm_mpz_t out, wm_mpz_t lhs, wm_mpz_t rhs);
void wm_mpz_divrem_raw_naive(wm_mpz_t div, wm_mpz_t rem, wm_mpz_t lhs, wm_mpz_t rhs);
void wm_mpz_divrem_raw_single_naive(wm_mpz_t div, wm_mpz_half_t* rem, wm_mpz_t lhs, wm_mpz_half_t rhs);

void wm_mpz_print_digits(wm_mpz_t z);
void wm_mpz_print(wm_mpz_t z);

#endif

#pragma once
#include <stdlib.h>
#include <malloc.h>
#include <intrin.h>
#include <string.h>

int silk_VAD_Get(
    //int          state,                       /*  Encoder state                               */
    const short            pIn[]                           /* I    PCM input                                   */
);

#define TYPE_NO_VOICE_ACTIVITY                  0
#define TYPE_UNVOICED                           1
#define TYPE_VOICED                             2

#define SPEECH_ACTIVITY_DTX_THRES                       0.05f
#define SILK_FIX_CONST( C, Q )              ((int)((C) * ((long)1 << (Q)) + 0.5))
#define silk_int16_MAX   0x7FFF                               /*  2^15 - 1 =  32767 */
#define silk_int16_MIN   ((short)0x8000)                 /* -2^15     = -32768 */
#define silk_int32_MAX   0x7FFFFFFF                           /*  2^31 - 1 =  2147483647 */
#define silk_int32_MIN   ((int)0x80000000)             /* -2^31     = -2147483648 */
#define silk_memset(dest, src, size)        memset((dest), (src), (size))

#define VAD_NOISE_LEVEL_SMOOTH_COEF_Q16         1024    /* Must be <  4096 */
#define VAD_NOISE_LEVELS_BIAS                   50

/* Sigmoid settings */
#define VAD_NEGATIVE_OFFSET_Q5                  128     /* sigmoid is 0 at -128 */
#define VAD_SNR_FACTOR_Q16                      45000

/* smoothing for SNR measurement */
#define VAD_SNR_SMOOTH_COEF_Q18                 4096

#define VAD_N_BANDS 4
#define VAD_INTERNAL_SUBFRAMES_LOG2             2
#define VAD_INTERNAL_SUBFRAMES                  ( 1 << VAD_INTERNAL_SUBFRAMES_LOG2 )
#define silk_uint8_MAX   0xFF                                 /*  2^8 - 1 = 255 */

#define VARDECL(type, var) type *var
#define silk_RSHIFT32(a, shift)             ((a)>>(shift))
#define silk_RSHIFT(a, shift)             ((a)>>(shift))
#define silk_LSHIFT32(a, shift)             ((a)<<(shift))
#define silk_LSHIFT(a, shift)             ((a)<<(shift))
#define ALLOC(var, size, type) var = ((type*)alloca(sizeof(type)*(size)))
#define silk_ADD16(a, b)                    ((a) + (b))
#define silk_ADD32(a, b)                    ((a) + (b))
#define silk_ADD64(a, b)                    ((a) + (b))

#define silk_SUB16(a, b)                    ((a) - (b))
#define silk_SUB32(a, b)                    ((a) - (b))
#define silk_SUB64(a, b)                    ((a) - (b))
#define silk_SMULWB(a32, b32)            ((((a32) >> 16) * (int)((short)(b32))) + ((((a32) & 0x0000FFFF) * (int)((short)(b32))) >> 16))
#define silk_SMLAWB(a32, b32, c32)       ((a32) + ((((b32) >> 16) * (int)((short)(c32))) + ((((b32) & 0x0000FFFF) * (int)((short)(c32))) >> 16)))
#define silk_SAT16(a)                       ((a) > silk_int16_MAX ? silk_int16_MAX :      \
                                            ((a) < silk_int16_MIN ? silk_int16_MIN : (a)))
#define silk_MLA(a32, b32, c32)             silk_ADD32((a32),((b32) * (c32)))
#define silk_SMLABB(a32, b32, c32)       ((a32) + ((int)((short)(b32))) * (int)((short)(c32)))
#define silk_ADD_POS_SAT32(a, b)            ((((unsigned int)(a)+(unsigned int)(b)) & 0x80000000) ? silk_int32_MAX : ((a)+(b)))
#define silk_ADD_POS_SAT32(a, b)            ((((unsigned int)(a)+(unsigned int)(b)) & 0x80000000) ? silk_int32_MAX : ((a)+(b)))
#define silk_DIV32_16(a32, b16)             ((int)((a32) / (b16)))
#define silk_DIV32(a32, b32)                ((int)((a32) / (b32)))
#define silk_RSHIFT_ROUND(a, shift)         ((shift) == 1 ? ((a) >> 1) + ((a) & 1) : (((a) >> ((shift) - 1)) + 1) >> 1)

#define silk_SMULWW(a32, b32)            silk_MLA(silk_SMULWB((a32), (b32)), (a32), silk_RSHIFT_ROUND((b32), 16))
#define silk_min(a, b)                      (((a) < (b)) ? (a) : (b))
#define silk_max(a, b)                      (((a) > (b)) ? (a) : (b))
#define silk_ADD_LSHIFT32(a, b, shift)      silk_ADD32((a), silk_LSHIFT32((b), (shift)))    /* shift >= 0 */
#define silk_MUL(a32, b32)                  ((a32) * (b32))
#define silk_SMULBB(a32, b32)            ((int)((short)(a32)) * (int)((short)(b32)))
#define silk_LIMIT( a, limit1, limit2)      ((limit1) > (limit2) ? ((a) > (limit1) ? (limit1) : ((a) < (limit2) ? (limit2) : (a))) \
                                                                 : ((a) > (limit2) ? (limit2) : ((a) < (limit1) ? (limit1) : (a))))

#define silk_LSHIFT_SAT32(a, shift)         (silk_LSHIFT32( silk_LIMIT( (a), silk_RSHIFT32( silk_int32_MIN, (shift) ), \
                                                    silk_RSHIFT32( silk_int32_MAX, (shift) ) ), (shift) ))







static const int tiltWeights[VAD_N_BANDS] = { 30000, 6000, -12000, -12000 };
static const int sigm_LUT_neg_Q15[6] = {
    16384, 8812, 3906, 1554, 589, 219
};
static const int sigm_LUT_slope_Q10[6] = {
    237, 153, 73, 30, 12, 7
};
static const int sigm_LUT_pos_Q15[6] = {
    16384, 23955, 28861, 31213, 32178, 32548
};

static __inline int ec_bsr(unsigned long _x) {
    unsigned long ret;
    _BitScanReverse(&ret, _x);
    return (int)ret;
}
# define EC_CLZ0    (1)
# define EC_CLZ(_x) (-ec_bsr(_x))
# define EC_ILOG(_x) (EC_CLZ0-EC_CLZ(_x))
static int silk_min_int(int a, int b)
{
    return (((a) < (b)) ? (a) : (b));
}
static int silk_max_int(int a, int b)
{
    return (((a) > (b)) ? (a) : (b));
}
static int silk_max_32(int a, int b)
{
    return (((a) > (b)) ? (a) : (b));
}
static  int silk_CLZ32(int in32)
{
    return in32 ? 32 - EC_ILOG(in32) : 32;
}
static  int silk_ROR32(int a32, int rot)
{
    unsigned int x = (unsigned int)a32;
    unsigned int r = (unsigned int)rot;
    unsigned int m = (unsigned int)-rot;
    if (rot == 0) {
        return a32;
    }
    else if (rot < 0) {
        return (int)((x << m) | (x >> (32 - m)));
    }
    else {
        return (int)((x << (32 - r)) | (x >> r));
    }
}
static  void silk_CLZ_FRAC(
    int in,            /* I  input                               */
    int* lz,           /* O  number of leading zeros             */
    int* frac_Q7       /* O  the 7 bits right after the leading one */
)
{
    int lzeros = silk_CLZ32(in);

    *lz = lzeros;
    *frac_Q7 = silk_ROR32(in, 24 - lzeros) & 0x7f;
}


/* Approximation of square root                                          */
/* Accuracy: < +/- 10%  for output values > 15                           */
/*           < +/- 2.5% for output values > 120                          */
static  int silk_SQRT_APPROX(int x)
{
    int y, lz, frac_Q7;

    if (x <= 0) {
        return 0;
    }

    silk_CLZ_FRAC(x, &lz, &frac_Q7);

    if (lz & 1) {
        y = 32768;
    }
    else {
        y = 46214;        /* 46214 = sqrt(2) * 32768 */
    }

    /* get scaling right */
    y >>= silk_RSHIFT(lz, 1);

    /* increment using fractional part of input */
    y = silk_SMLAWB(y, y, silk_SMULBB(213, frac_Q7));

    return y;
}
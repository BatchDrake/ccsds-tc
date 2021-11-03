/*

  Copyright (C) 2021 Gonzalo José Carracedo Carballal

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, version 3.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this program.  If not, see
  <http://www.gnu.org/licenses/>

*/

#ifndef _CCSDS_TC_H
#define _CCSDS_TC_H

#define _ISOC99_SOURCE

#include <stdbool.h>
#include <math.h>
#include <stdint.h>

#define CCSDS_TC_MINIMUM_BLOCK_LENGTH 223

#ifdef __GNUC__
#  define CCSDS_TC_INLINE static inline __attribute__((always_inline))
#else
#  define CCSDS_TC_INLINE static inline
#endif

#ifdef CCSDS_TC_INT_ARITHMETICS
#  define CCSDS_TC_SQRT_BITS 6
#  define CCSDS_TC_SSQRT (1 << CCSDS_TC_SQRT_BITS)
#  define CCSDS_TC_SCALE (1 << (CCSDS_TC_SQRT_BITS << 1))
#  define CCSDS_TC_MAX_LH  0
#  define CCSDS_TC_MIN_LH  ((softbit_t) 0x80000000)
#  define CCSDS_TC_IMPOSSIBLE(x) ((x) == 0x80000000)
#  define CCSDS_TC_MAX(a, b) ((a) > (b) ? (a) : (b))
#  define CCSDS_TC_SAFE_SUM(a, b)                   \
  ((a) == CCSDS_TC_MIN_LH || (b) == CCSDS_TC_MIN_LH \
    ? CCSDS_TC_MIN_LH                               \
    : (a) + (b))

typedef int32_t softbit_t;
#else
#  define CCSDS_TC_SCALE 1.0
#  define CCSDS_TC_SSQRT 1.0
#  define CCSDS_TC_MAX_LH  0
#  define CCSDS_TC_MIN_LH -INFINITY
#  define CCSDS_TC_IMPOSSIBLE(x) isinf(x)
#  define CCSDS_TC_MAX(a, b) fmax(a, b)
#  define CCSDS_TC_SAFE_SUM(a, b) ((a) + (b))
typedef float softbit_t;
#endif


#define CCSDS_TC_F2SB(x) ((softbit_t) ((x) * CCSDS_TC_SCALE))
#define CCSDS_TC_SB2F(x) ((float) (x) / (float) CCSDS_TC_SCALE)

#define CCSDS_TC_STATE_BITS 4

#define CCSDS_TC_STATE_NUM    (1 << CCSDS_TC_STATE_BITS)
#define CCSDS_TC_STATE_MASK   (CCSDS_TC_STATE_NUM - 1)
#define CCSDS_TC_STATE_ZERO   0
#define CCSDS_TC_STATE_UNUSED 0xff

#define CCSDS_TC_G0 0x19 /* 11001 */
#define CCSDS_TC_G1 0x1b /* 11011 */
#define CCSDS_TC_G2 0x15 /* 10101 */
#define CCSDS_TC_G3 0x1f /* 11111 */

#define CCSDS_TC_EBN0_TO_LC(ebn0) (4 * (ebn0))

/* Better when state matches word length */
typedef uint32_t state_t;
typedef state_t  poly_t;

/************************* CONVOLUTIONAL CODE API ***************************/

/*
 *  HOW IS THE CONVOLUTIONAL CODE ACTUALLY REPRESENTED
 *
 * ccsds_tc_bcjr assumes a generic convolutional code
 * over a recursive shift register with the following
 * internal wiring:
 * 
 *     cn-1      c3      c2      c1       c0    
 *  Dn-1 <- ... <--- D2 <--- D1 <--- D0 <---(+)---- u 
 *  | cn             |       |       |       ^
 *  v                v       v       v       |
 * (?)----- ... ----(?)-----(?)-----(?)------/
 * bn               b3      b2      b1
 * 
 * In order to recover a systematic output, the output
 * of the XOR right before c0 must be XOR'ed back by the
 * backward connections again. This is easily achiveable
 * if the systematic forward output polynomial matches
 * the backward polynomial up to c1. c0 must be set to 1
 * in this case.
 *  
 * On the other hand, the j-th bit of the backward polynomial
 * refers to the output of the j-1-th latch. b1 is connected to
 * the output of D0, b2 is connected to the otput of D1, etc.
 * 
 * For instance, the example provided by Ammar Abh which was based
 * on a systematic convolutional code of the following structure:
 * 
 *     c0     
 * D0 <---(+)--- u
 * |       ^          v0 = c0 ^ c1
 * | c1    |          v1 = c0
 * \-------/
 * 
 * The forward polynomials for v0 and v1 are:
 *   G(v0) = 11
 *   G(v1) = 01
 *   
 * While the backward connection polynomial is simply:
 *   G(bi) = 1X
 * 
 * The state_t type represents the ltches D0...Dn-1
 */

struct ccsds_tc_bcjr_ctx {
  unsigned int h;      /* Information bits */
  unsigned int length; /* Data size as accepted by the algo */
  unsigned int n;      /* Codeword length */

  /* Cached data */
  softbit_t *rv;
  bool rv_cached;
  
  state_t S_last_0;
  state_t S_last_1;
  
  /* Trellis connections */
  poly_t     poly_b;
  state_t   *S_next;
  state_t   *S_prev;

  /* Codeword table */
  state_t   *hv;

  /* Path metrics */
  softbit_t *astar;
  softbit_t *bstar;

  /* Branch metric */

  /*
   * gstar_dir[2 * (STATE_NUM * l + Si) + 0]: gamma_l*(Si, next(Si, 0))
   * gstar_dir[2 * (STATE_NUM * l + Si) + 1]: gamma_l*(Si, next(Si, 1))
   */
  softbit_t *gstar_dir;

  /*
   * gstar_rev[2 * (STATE_NUM * l + Si) + 0]: gamma_l*(prev(Si, 0), Si)
   * gstar_rev[2 * (STATE_NUM * l + Si) + 1]: gamma_l*(prev(Si, 1), Si)
   */
  softbit_t *gstar_rev;
  softbit_t *min_lh_buf;
};

typedef struct ccsds_tc_bcjr_ctx ccsds_tc_bcjr_ctx_t;

bool ccsds_tc_bcjr_ctx_init(
  ccsds_tc_bcjr_ctx_t *self,
  const poly_t *poly_f,
  poly_t poly_b,
  unsigned int n,
  unsigned int h);

bool ccsds_tc_bcjr_ctx_decode(
  ccsds_tc_bcjr_ctx_t *self,
  const softbit_t *x,
  const softbit_t *L_u,
  softbit_t *y,
  softbit_t L_c,
  const int *restrict permute_out);

CCSDS_TC_INLINE void
ccsds_tc_bcjr_invalidate_cache(ccsds_tc_bcjr_ctx_t *self)
{
  self->rv_cached = false;
}

void ccsds_tc_bcjr_ctx_finalize(ccsds_tc_bcjr_ctx_t *self);

/************************* CONVOLUTIONAL CODE API ***************************/
enum ccsds_tc_decoder_block_length {
  CCSDS_TC_BLOCK_LENGTH_223  = 1,
  CCSDS_TC_BLOCK_LENGTH_446  = 2,
  CCSDS_TC_BLOCK_LENGTH_892  = 4,
  CCSDS_TC_BLOCK_LENGTH_1115 = 5
};

struct ccsds_tc_decoder_params {
  unsigned int iters;
  unsigned int rate;
  unsigned int block_len;
  softbit_t    L_c;
};

#define ccsds_tc_decoder_params_INITIALIZER  \
{                                            \
  20, /* iters */                            \
  6,  /* rate */                             \
  CCSDS_TC_BLOCK_LENGTH_223, /* block_len */ \
  4, /* L_c */                               \
}

struct ccsds_tc_decoder {
  struct ccsds_tc_decoder_params params;
  softbit_t *x1, *x2;
  softbit_t *y1, *y2;
  softbit_t *y0;

  void (*ser2par_sb_cb) (struct ccsds_tc_decoder *, const softbit_t *);
  void (*ser2par_float_cb) (struct ccsds_tc_decoder *, const float *);

  int *permutate_dir;
  int *permutate_inv;
  unsigned int n1;
  unsigned int rate;
  unsigned int length;
  unsigned int terminated_len;
  unsigned int cw_length;
  ccsds_tc_bcjr_ctx_t dec1;
  ccsds_tc_bcjr_ctx_t dec2;
};

typedef struct ccsds_tc_decoder ccsds_tc_decoder_t;

bool ccsds_tc_decoder_init(
  ccsds_tc_decoder_t *self, 
  const struct ccsds_tc_decoder_params *params);

CCSDS_TC_INLINE unsigned int
ccsds_tc_decoder_get_block_size(const ccsds_tc_decoder_t *self)
{
  return self->length;
}

CCSDS_TC_INLINE unsigned int
ccsds_tc_decoder_get_codeword_size(const ccsds_tc_decoder_t *self)
{
  return self->cw_length;
}

CCSDS_TC_INLINE float
ccsds_tc_get_error_rate(const ccsds_tc_decoder_t *self)
{
  unsigned int i, count = 0;

  for (i = 0; i < self->length; ++i)
    count += (self->x1[self->n1 * i] < 0) == (self->y2[i] < 0);

  return (float) count / (float) self->length;
}

CCSDS_TC_INLINE const softbit_t *
ccsds_tc_get_intermediate_output(const ccsds_tc_decoder_t *self)
{
  return self->y1;
}

CCSDS_TC_INLINE const softbit_t *
ccsds_tc_get_output(const ccsds_tc_decoder_t *self)
{
  return self->y2;
}

void ccsds_tc_decoder_feed_block(ccsds_tc_decoder_t *self, const softbit_t *block);
void ccsds_tc_decoder_feed_block_float(ccsds_tc_decoder_t *self, const float *block);

void ccsds_tc_decoder_finalize(ccsds_tc_decoder_t *self);

#endif /* _CCSDS_TC_H */

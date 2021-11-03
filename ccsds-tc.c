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

#include "ccsds-tc.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define LOG2              CCSDS_TC_F2SB(0.69314718055995)
#define MAXSTAR_LIN_LIMIT CCSDS_TC_F2SB(2 * 0.69314718055995)
#define BIT2SOFT(b) (CCSDS_TC_SCALE * (2 * (softbit_t) (b) - 1))
#define ALLOC(where, N, what)                     \
  if ((where = malloc(N * sizeof(what))) == NULL) \
    goto fail


/******************* BCJR convolutional decoder *********************/
#ifdef CCSDS_TC_NO_MAXSTAR
#  define ccsds_tc_softbit_maxstar CCSDS_TC_MAX
#else
CCSDS_TC_INLINE softbit_t
ccsds_tc_softbit_maxstar(softbit_t a, softbit_t b)
{
  softbit_t diff;

  if (CCSDS_TC_IMPOSSIBLE(a) || CCSDS_TC_IMPOSSIBLE(b))
    return CCSDS_TC_MAX(a, b);

  if (a > b) {
    diff = a - b;
    if (diff > MAXSTAR_LIN_LIMIT)
      return a;
    else
      return a + LOG2 - diff / 2;
  } else {
    diff = b - a;
    if (diff > MAXSTAR_LIN_LIMIT)
      return b;
    else
      return b + LOG2 - diff / 2;
  }
}
#endif /* CCSDS_TC_NO_MAXSTAR */

/* This is derived from a Taylor expansion of
     fmaxf(a, b) + log(1 + exp(-fabs(a-b))) */

CCSDS_TC_INLINE state_t
ccsds_tc_bcjr_state_dotprod(state_t a, state_t b) 
{
#ifdef __GNUC__
  return __builtin_parity(a & b);
#else
  /* TODO: generate at compile time */
  return
      ((a >> 0) & (b >> 0) & 1)
    ^ ((a >> 1) & (b >> 1) & 1)
    ^ ((a >> 2) & (b >> 2) & 1)
    ^ ((a >> 3) & (b >> 3) & 1)
    ^ ((a >> 4) & (b >> 4) & 1);
#endif  
}

CCSDS_TC_INLINE state_t
ccsds_tc_bcjr_state_next(state_t a, state_t poly_b, bool next)
{
  next ^= ccsds_tc_bcjr_state_dotprod(a, poly_b);

  return CCSDS_TC_STATE_MASK & ((a << 1) | (next & 1));
}

CCSDS_TC_INLINE bool
ccsds_tc_bcjr_connect_state(
  ccsds_tc_bcjr_ctx_t *self,
  state_t prev,
  state_t next,
  bool b)
{
  /* Forward connection. Trivial. */
  self->S_next[(prev << 1) + b] = next;

  /* 
   * Backward connection. Not so trivial. We check whether
   * the previous state is set to "unused" and update the
   * entry accordingly
   */

  if (next == 0) {
    if (b)
      self->S_last_1 = prev;
    else
      self->S_last_0 = prev;
  }

  if (self->S_prev[(next << 1) + 0] == CCSDS_TC_STATE_UNUSED) {
    self->S_prev[(next << 1) + 0] = prev;
  } else if (self->S_prev[(next << 1) + 1] == CCSDS_TC_STATE_UNUSED) {
    self->S_prev[(next << 1) + 1] = prev;
  } else {
    fprintf(stderr, "critical: more than 2 previous states?\n");
    return false;
  }

  return true;
}

void
ccsds_tc_bcjr_ctx_finalize(ccsds_tc_bcjr_ctx_t *self)
{
  if (self->S_next != NULL)
    free(self->S_next);

  if (self->S_prev != NULL)
    free(self->S_prev);

  if (self->hv != NULL)
    free(self->hv);

  if (self->rv != NULL)
    free(self->rv);

  if (self->astar != NULL)
    free(self->astar);

  if (self->bstar != NULL)
    free(self->bstar);

  if (self->gstar_dir != NULL)
    free(self->gstar_dir);

  if (self->gstar_rev != NULL)
    free(self->gstar_rev);

  if (self->min_lh_buf != NULL)
    free(self->min_lh_buf);
}

bool
ccsds_tc_bcjr_ctx_init(
  ccsds_tc_bcjr_ctx_t *self,
  const poly_t *poly_f,
  poly_t poly_b,
  unsigned int n,
  unsigned int h)
{
  unsigned int i, j;
  unsigned int cwnb; /* Codeword number and bit */
  unsigned int cwnb_0, cwnb_1;
  bool b;
  state_t next0, next1;
  softbit_t lh;
  bool ok = false;

  memset(self, 0, sizeof(ccsds_tc_bcjr_ctx_t));

  poly_b >>= 1; /* Discard rightmost bit */

  self->h = h;
  self->n = n;
  self->poly_b = poly_b;
  self->length = h + CCSDS_TC_STATE_BITS ;

  /* Allocate data */
  ALLOC(self->hv,               2 * CCSDS_TC_STATE_NUM,      state_t);
  ALLOC(self->S_next, 2 * self->n * CCSDS_TC_STATE_NUM,      state_t);
  ALLOC(self->S_prev, 2 * self->n * CCSDS_TC_STATE_NUM,      state_t);
  ALLOC(self->rv,     2 * self->length * CCSDS_TC_STATE_NUM, softbit_t);
  ALLOC(self->astar,      self->length * CCSDS_TC_STATE_NUM, softbit_t);
  ALLOC(self->bstar,      (self->length + 1) * CCSDS_TC_STATE_NUM, softbit_t);
  ALLOC(self->gstar_dir,  2 * self->length * CCSDS_TC_STATE_NUM,   softbit_t);
  ALLOC(self->gstar_rev,  2 * self->length * CCSDS_TC_STATE_NUM,   softbit_t);
  ALLOC(self->min_lh_buf, 2 * CCSDS_TC_STATE_NUM,                  softbit_t);

  for (j = 0; j < 2 * self->length * CCSDS_TC_STATE_NUM; ++j)
    self->gstar_dir[j] = self->gstar_rev[j] = CCSDS_TC_MIN_LH;
  
  for (j = 0; j < 2 * CCSDS_TC_STATE_NUM; ++j)
    self->min_lh_buf[j] = CCSDS_TC_MIN_LH;
  
  memset(self->hv, 0, 2 * CCSDS_TC_STATE_NUM * sizeof(state_t));
  memset(self->rv, 0, 2 * self->length * CCSDS_TC_STATE_NUM * sizeof(softbit_t));

  /* Mark all previous states as unused */
  for (j = 0; j < CCSDS_TC_STATE_NUM; ++j) {
    self->S_prev[(j << 1) + 0] = CCSDS_TC_STATE_UNUSED;
    self->S_prev[(j << 1) + 1] = CCSDS_TC_STATE_UNUSED;
  }
  
  /* Initialize path metrics */
  for (j = 0; j < CCSDS_TC_STATE_NUM; ++j) {
    lh = j == 0 ? CCSDS_TC_MAX_LH : CCSDS_TC_MIN_LH;
    self->astar[j] = lh;
    self->bstar[CCSDS_TC_STATE_NUM * (self->length) + j] = lh;
  }

  /* Initialize codeword table AND trellis connections. -1 is zero. */
  for (j = 0; j < CCSDS_TC_STATE_NUM; ++j) {
    b     = ccsds_tc_bcjr_state_dotprod(j, poly_b);
    next0 = ccsds_tc_bcjr_state_next(j, poly_b, false);
    next1 = ccsds_tc_bcjr_state_next(j, poly_b, true);

    /* Connect state j to next if input is 0 */
    if (!ccsds_tc_bcjr_connect_state(
      self,
      j,
      next0,
      false))
      goto fail;

    /* Connect state j to next if input is 1 */
    if (!ccsds_tc_bcjr_connect_state(
      self,
      j,
      next1,
      true))
      goto fail;

    for (i = 0; i < n; ++i)
      self->hv[(j << 1) + 0] 
        |= ccsds_tc_bcjr_state_dotprod(poly_f[i], (j << 1) ^ b ^ false) << i;
    for (i = 0; i < n; ++i)
      self->hv[(j << 1) + 1] 
        |= ccsds_tc_bcjr_state_dotprod(poly_f[i], (j << 1) ^ b ^ true) << i;
  }

  ok = true;

fail:  
  if (!ok)
    ccsds_tc_bcjr_ctx_finalize(self);

  return ok;
}

CCSDS_TC_INLINE softbit_t
ccsds_tc_bcjr_ctx_semisoft_dotprod_2(const softbit_t *__restrict r, state_t hv)
{
  switch (hv) {
    case 0:
      return -(r[0] + r[1]);
    case 1:
      return r[0] - r[1];
    case 2:
      return r[1] - r[0];
    case 3:
      return r[1] + r[0];
  }

  return 0;
}
CCSDS_TC_INLINE softbit_t
ccsds_tc_bcjr_ctx_semisoft_dotprod(const softbit_t *__restrict r, state_t hv, unsigned int L)
{
  softbit_t accum = 0;

  switch (L) {
    case 2:
      return ccsds_tc_bcjr_ctx_semisoft_dotprod_2(r, hv);

    case 4:
      return ccsds_tc_bcjr_ctx_semisoft_dotprod_2(r, hv & 3)
      + ccsds_tc_bcjr_ctx_semisoft_dotprod_2(r + 2, (hv >> 2) & 3);
  }

  while (L-- != 0) {
    if (hv & 1)
      accum += *r++;
    else
      accum -= *r++;

    hv >>= 1;
  }
  return accum;
}


CCSDS_TC_INLINE void
ccsds_tc_bcjr_ctx_set_gamma(
  ccsds_tc_bcjr_ctx_t *self,
  unsigned int j, /* position */
  state_t state,  /* previous state */
  bool b,         /* next bit */
  const softbit_t *__restrict L_u)   /* received codeword */
{
  unsigned int sndx = (state << 1) + b;
  unsigned int pos_off = j * CCSDS_TC_STATE_NUM << 1;
  unsigned int dest = pos_off + sndx;
  state_t next = self->S_next[sndx];
  state_t next2 = next << 1;

 /*
  * Gamma star is already a posterior log-probability
  * that the received codeword corresponds to a given bit
  * given the current state. It is an expression of the form:
  * 
  * 1/2 (bit * L(bit) + CW * L(CW))
  * 
  * In this expression we basically combine prior information
  * of the bit (L(bit)) + the soft information of the received
  * codeword. This is dot-multiplied by the bit and expected
  * codeword of the current branch.
  */

  /* We only set gamma in the intermediate cases */
  const softbit_t gstar = (j < self->h 
    ? (BIT2SOFT(b) / CCSDS_TC_SSQRT) * (L_u[j] / (2 * CCSDS_TC_SSQRT)) : 0) 
    + self->rv[dest];

  /* 
  * For each j, gamma is actually a multivaluated function
  * on the state number. Since we don't actually need
  * STATE_NUM x STATE_NUM gammas, as many of the state connections
  * are forbidden, we provide a sparse representation in
  * which the row-to-cols and col-to-rows view are stored.
  * We can do this because we know that, for each departure
  * state there are 2 possible arrival states (gamma_dir) and
  * for each arrival states, there are 2 possible departure
  * states (gamma_rev).
  */

  /* gstar_dir*: gamma* row-to-cols view */
  self->gstar_dir[dest] = gstar;
  
  /* 
  * gstar_rev: gamma* col-to-rows view. Now this is
  * the trick: we are going to enforce the SAME indexing
  * as with S_prev, so we know that the first of the
  * two gamma* arriving to state "next" refers to the
  * the first of the previous states to "next" 
  */
  
  self->gstar_rev[pos_off + next2 + (self->S_prev[next2] != state)] = gstar;
}


CCSDS_TC_INLINE void
ccsds_tc_bcjr_ctx_set_gamma_bulk(
  ccsds_tc_bcjr_ctx_t *self,
  const softbit_t *__restrict L_u)
{
  unsigned int p = CCSDS_TC_STATE_NUM << 1;
  unsigned int i, j;

#pragma GCC unroll 2
  for (j = 1; j < self->h; ++j) {
#pragma GCC unroll 32
    for (i = 0; i < (CCSDS_TC_STATE_NUM << 1); ++i, ++p) {
      const state_t next    = self->S_next[i];
      const softbit_t gstar = 
        (BIT2SOFT(i & 1) / CCSDS_TC_SSQRT) * (L_u[j] / (2 * CCSDS_TC_SSQRT))
        + self->rv[p];
      const unsigned int rev_p = 
        ((j * CCSDS_TC_STATE_NUM + next) << 1) + (self->S_prev[next << 1] != (i >> 1));
      self->gstar_dir[p]     = gstar;
      self->gstar_rev[rev_p] = gstar;
    }
  }

#pragma GCC unroll 2
  for (j = self->h; j < self->length - 1; ++j) {
#pragma GCC unroll 32
    for (i = 0; i < (CCSDS_TC_STATE_NUM << 1); ++i, ++p) {
      const state_t next    = self->S_next[i];
      const softbit_t gstar = self->rv[p];
      const unsigned int rev_p = 
        ((j * CCSDS_TC_STATE_NUM + next) << 1) + (self->S_prev[next << 1] != (i >> 1));
      self->gstar_dir[p]     = gstar;
      self->gstar_rev[rev_p] = gstar;
    }
  }
}

CCSDS_TC_INLINE void
ccsds_tc_bcjr_ctx_set_alpha(
  ccsds_tc_bcjr_ctx_t *self,
  unsigned int j, /* next position */
  state_t next)   /* next state */
{
  unsigned int pos_off  = j * CCSDS_TC_STATE_NUM;
  unsigned int prev_off = pos_off - CCSDS_TC_STATE_NUM;
  unsigned int gprv_off = prev_off << 1;
  unsigned int i0 = (next << 1) + 0;
  unsigned int i1 = i0 + 1;
  softbit_t astar;

  /* There are only two previous states */
  astar = 
    ccsds_tc_softbit_maxstar(
      CCSDS_TC_SAFE_SUM(
        self->gstar_rev[gprv_off + i0], 
        self->astar[prev_off + self->S_prev[i0]]),
      CCSDS_TC_SAFE_SUM(
        self->gstar_rev[gprv_off + i1], 
        self->astar[prev_off + self->S_prev[i1]]));

  self->astar[pos_off + next] = astar;
}

CCSDS_TC_INLINE void
ccsds_tc_bcjr_ctx_set_beta(
  ccsds_tc_bcjr_ctx_t *self,
  unsigned int j, /* curr position */
  state_t prev)   /* curr state */
{
  unsigned int pos_off = j * CCSDS_TC_STATE_NUM;
  unsigned int next_off = pos_off + CCSDS_TC_STATE_NUM;
  unsigned int gpos_off = pos_off << 1;
  unsigned int i0 = (prev << 1) + 0;
  unsigned int i1 = i0 + 1;
  softbit_t bstar;

  /* There are only two next states */
  bstar = 
    ccsds_tc_softbit_maxstar(
      CCSDS_TC_SAFE_SUM(
        self->gstar_dir[gpos_off + i0],
        self->bstar[next_off + self->S_next[i0]]),
      CCSDS_TC_SAFE_SUM(
        self->gstar_dir[gpos_off + i1],
        self->bstar[next_off + self->S_next[i1]]));

  self->bstar[pos_off + prev] = bstar;
}

CCSDS_TC_INLINE softbit_t
ccsds_tc_bcjr_ctx_get_Lu(ccsds_tc_bcjr_ctx_t *self, unsigned int j)
{
  unsigned int pos_off  = j * CCSDS_TC_STATE_NUM;
  unsigned int pos_off2 = pos_off << 1;
  unsigned int next_off = pos_off + CCSDS_TC_STATE_NUM;
  unsigned int i, state_off;
  softbit_t L0, L1;

  L0 = CCSDS_TC_MIN_LH;
  L1 = CCSDS_TC_MIN_LH;

  state_off = 0;
#pragma GCC unroll 16
  for (i = 0; i < CCSDS_TC_STATE_NUM; ++i) {  
    L0 = ccsds_tc_softbit_maxstar(
      CCSDS_TC_SAFE_SUM(
        self->bstar[next_off + self->S_next[state_off]],
        CCSDS_TC_SAFE_SUM(
          self->gstar_dir[pos_off2 + state_off],
        self->astar[pos_off + i])),
      L0);
    state_off += 2;
  }
  
  state_off = 1;

#pragma GCC unroll 16
  for (i = 0; i < CCSDS_TC_STATE_NUM; ++i) {
    L1 = ccsds_tc_softbit_maxstar(
      CCSDS_TC_SAFE_SUM(
        self->bstar[next_off + self->S_next[state_off]],
        CCSDS_TC_SAFE_SUM(
          self->gstar_dir[pos_off2 + state_off],
          self->astar[pos_off + i])),
      L1);
    state_off += 2;
  }

  return L1 - L0;
}

CCSDS_TC_INLINE softbit_t
ccsds_tc_bcjr_ctx_rebuild_rv_cache(
  ccsds_tc_bcjr_ctx_t *self,
  const softbit_t *__restrict x,
  softbit_t L_c)
{
  unsigned int i, j, p = 0;
  unsigned int q = 0;

#pragma GCC unroll 2
  for (j = 0; j < self->length; ++j) {
#pragma GCC unroll 32
    for (i = 0; i < (CCSDS_TC_STATE_NUM << 1); ++i, ++p) {
      self->rv[p] 
        = (L_c / (2 * CCSDS_TC_SSQRT))
        * (ccsds_tc_bcjr_ctx_semisoft_dotprod(
          &x[q], 
          self->hv[i], 
          self->n) / CCSDS_TC_SSQRT);
    }
    q += self->n;
  }

  self->rv_cached = true;
}

bool
ccsds_tc_bcjr_ctx_decode(
  ccsds_tc_bcjr_ctx_t *self,
  const softbit_t *__restrict x,
  const softbit_t *__restrict L_u,
  softbit_t *y,
  softbit_t L_c,
  const int *__restrict permute_out)
{
  unsigned int i, j, i0, i1, off = 0;
  softbit_t lh, L0, L1;
  unsigned int cwnb;
  unsigned int cwnb_0, cwnb_1;
  bool ok = false;

  if (!self->rv_cached)
    ccsds_tc_bcjr_ctx_rebuild_rv_cache(self, x, L_c);

  ccsds_tc_bcjr_ctx_set_gamma(self, 0, 0, false, L_u);
  ccsds_tc_bcjr_ctx_set_gamma(self, 0, 0, true,  L_u);
  ccsds_tc_bcjr_ctx_set_gamma(
    self, 
    self->length - 1,
    self->S_last_0,  
    false, 
    L_u);
  ccsds_tc_bcjr_ctx_set_gamma(
    self, 
    self->length - 1,
    self->S_last_1,  
    true, 
    L_u);

  ccsds_tc_bcjr_ctx_set_gamma_bulk(self, L_u);

  /* Set path metrics (forward) */
#pragma GCC unroll 2
  for (j = 1; j < self->length; ++j)
#pragma GCC unroll 16
    for (i = 0; i < CCSDS_TC_STATE_NUM; ++i)
      ccsds_tc_bcjr_ctx_set_alpha(self, j, i);

  /* Set path metrics (backward) */
#pragma GCC unroll 2
  for (j = 1; j < self->length; ++j)
#pragma GCC unroll 16
    for (i = 0; i < CCSDS_TC_STATE_NUM; ++i)
      ccsds_tc_bcjr_ctx_set_beta(self, self->length - j, i);

  /* Compute bit likelihoods and finish */
  if (permute_out != NULL) {
#pragma GCC unroll 2
  for (j = 0; j < self->h; ++j)
    y[j] = ccsds_tc_bcjr_ctx_get_Lu(self, permute_out[j]);
#pragma GCC unroll 22
  for (j = self->h; j < self->length; ++j)
    y[j] = 0;
  } else {
#pragma GCC unroll 2
  for (j = 0; j < self->length; ++j)
    y[j] = ccsds_tc_bcjr_ctx_get_Lu(self, j);
  }

  ok = true;

done:
  return ok;
}

/*********************** Full CCSDS Turbocode *******************/
static int *
ccsds_tc_decoder_make_permutator(int k2)
{
  int *pi = NULL;
  int s, m, i, j, t, q, c;
  int p[] = {31, 37, 43, 47, 53, 59, 61, 67};

  if ((pi = malloc(sizeof(int) * k2 * 8)) == NULL)
    return NULL;

  for (s = 1; s <= k2 * 8; ++s) {
    m = (s - 1) & 1;
    i = (s - 1) / (2 * k2);
    j = (s - 1) / 2 - i * k2;
    t = (19 * i + 1) % 4;
    q = (t % 8) + 1;
    c = (p[q - 1] * j + 21 * m) % k2;
    pi[s - 1] = 2 * (t + c * 4 + 1) - m - 1;
  }

  return pi;
}

/* RATE 1/2: Handmade puncturing */
static void
ser2par_sb_2(ccsds_tc_decoder_t *self, const softbit_t *__restrict r)
{
  unsigned int i;
  
  for (i = 0; i < self->terminated_len; ++i) {
    self->x1[2 * i + 0] = r[2 * i + 0];
    self->x1[2 * i + 1] = i & 1 ? CCSDS_TC_F2SB(0) : r[2 * i + 1];
    self->x2[i]         = i & 1 ? r[2 * i + 1]     : CCSDS_TC_F2SB(0);
  }
}

static void
ser2par_float_2(ccsds_tc_decoder_t *self, const float *__restrict r)
{
  unsigned int i;
  
  for (i = 0; i < self->terminated_len; ++i) {
    self->x1[2 * i + 0] = CCSDS_TC_F2SB(r[2 * i + 0]);
    self->x1[2 * i + 1] = i & 1 ? CCSDS_TC_F2SB(0) : CCSDS_TC_F2SB(r[2 * i + 1]);
    self->x2[i] = i & 1 ? CCSDS_TC_F2SB(r[2 * i + 1]) : CCSDS_TC_F2SB(0);
  }
}

/* RATE 1/3 */
static void
ser2par_sb_3(ccsds_tc_decoder_t *self, const softbit_t *__restrict r)
{
  unsigned int i;
  
  for (i = 0; i < self->terminated_len; ++i) {
    self->x1[2 * i + 0] = r[3 * i + 0];
    self->x1[2 * i + 1] = r[3 * i + 1];
    self->x2[i]         = r[3 * i + 2];
  }
}

static void
ser2par_float_3(ccsds_tc_decoder_t *self, const float *__restrict r)
{
  unsigned int i;
  
  for (i = 0; i < self->terminated_len; ++i) {
    self->x1[2 * i + 0] = CCSDS_TC_F2SB(r[3 * i + 0]);
    self->x1[2 * i + 1] = CCSDS_TC_F2SB(r[3 * i + 1]);
    self->x2[i]         = CCSDS_TC_F2SB(r[3 * i + 2]);
  }
}

/* RATE 1/4 */
static void
ser2par_sb_4(ccsds_tc_decoder_t *self, const softbit_t *__restrict r)
{
  unsigned int i;
  
  for (i = 0; i < self->terminated_len; ++i) {
    self->x1[3 * i + 0] = r[4 * i + 0];
    self->x1[3 * i + 1] = r[4 * i + 1];
    self->x1[3 * i + 2] = r[4 * i + 2];
    self->x2[i]         = r[4 * i + 3];
  }
}

static void
ser2par_float_4(ccsds_tc_decoder_t *self, const float *__restrict r)
{
  unsigned int i;
  
  for (i = 0; i < self->terminated_len; ++i) {
    self->x1[3 * i + 0] = CCSDS_TC_F2SB(r[4 * i + 0]);
    self->x1[3 * i + 1] = CCSDS_TC_F2SB(r[4 * i + 1]);
    self->x1[3 * i + 2] = CCSDS_TC_F2SB(r[4 * i + 2]);
    self->x2[i]         = CCSDS_TC_F2SB(r[4 * i + 3]);
  }
}

/* RATE 1/6 */
static void
ser2par_sb_6(ccsds_tc_decoder_t *self, const softbit_t *__restrict r)
{
  unsigned int i;
  
  for (i = 0; i < self->terminated_len; ++i) {
    self->x1[4 * i + 0] = r[6 * i + 0];
    self->x1[4 * i + 1] = r[6 * i + 1];
    self->x1[4 * i + 2] = r[6 * i + 2];
    self->x1[4 * i + 3] = r[6 * i + 3];
    self->x2[2 * i + 0] = r[6 * i + 4];
    self->x2[2 * i + 1] = r[6 * i + 5];
  }
}

static void
ser2par_float_6(ccsds_tc_decoder_t *self, const float *__restrict r)
{
  unsigned int i;
  
  for (i = 0; i < self->terminated_len; ++i) {
    self->x1[4 * i + 0] = CCSDS_TC_F2SB(r[6 * i + 0]);
    self->x1[4 * i + 1] = CCSDS_TC_F2SB(r[6 * i + 1]);
    self->x1[4 * i + 2] = CCSDS_TC_F2SB(r[6 * i + 2]);
    self->x1[4 * i + 3] = CCSDS_TC_F2SB(r[6 * i + 3]);
    self->x2[2 * i + 0] = CCSDS_TC_F2SB(r[6 * i + 4]);
    self->x2[2 * i + 1] = CCSDS_TC_F2SB(r[6 * i + 5]);
  }
}

bool
ccsds_tc_decoder_init(
  ccsds_tc_decoder_t *self, 
  const struct ccsds_tc_decoder_params *params)
{
  unsigned int block_len;
  unsigned int block_bits;
  unsigned int i;
  unsigned int terminated_len;
  unsigned int n1[] = {0, 2, 2, 3, 0, 4};
  unsigned int n2[] = {0, 1, 1, 1, 0, 2};
  poly_t poly_b  = CCSDS_TC_G0;
  poly_t poly_f1[4];
  poly_t poly_f2[2];
  
  unsigned int rate;

  bool ok = false;

  memset(self, 0, sizeof(ccsds_tc_decoder_t));

  self->params = *params;

  /* Check block length */
  if (self->params.block_len < CCSDS_TC_BLOCK_LENGTH_223 
  || self->params.block_len > CCSDS_TC_BLOCK_LENGTH_1115
  || self->params.block_len == 3) {
    fprintf(stderr, "%s: invalid block length\n", __FUNCTION__);
    goto fail;
  }

  /* Configure channel rate */
  poly_f1[0] = CCSDS_TC_G0; /* The same for all rates */
  poly_f2[0] = CCSDS_TC_G1; /* Ditto */

  switch (self->params.rate) {
    case 2:
      self->ser2par_sb_cb    = ser2par_sb_2;
      self->ser2par_float_cb = ser2par_float_2;
      poly_f1[1]             = CCSDS_TC_G1;
      break;
    
    case 3:
      self->ser2par_sb_cb    = ser2par_sb_3;
      self->ser2par_float_cb = ser2par_float_3;
      poly_f1[1]             = CCSDS_TC_G1;
      break;
    
    case 4:
      self->ser2par_sb_cb    = ser2par_sb_4;
      self->ser2par_float_cb = ser2par_float_4;
      poly_f1[1]             = CCSDS_TC_G2;
      poly_f1[2]             = CCSDS_TC_G3;
      break;
    
    case 6:
      self->ser2par_sb_cb    = ser2par_sb_6;
      self->ser2par_float_cb = ser2par_float_6;
      poly_f1[1]             = CCSDS_TC_G1;
      poly_f1[2]             = CCSDS_TC_G2;
      poly_f1[3]             = CCSDS_TC_G3;
      poly_f2[1]             = CCSDS_TC_G3;
      break;
    
    default:
      fprintf(stderr, "%s: invalid rate 1/%d\n", self->params.rate);
      goto fail;
  }

  rate           = self->params.rate;
  block_len      = self->params.block_len * CCSDS_TC_MINIMUM_BLOCK_LENGTH;
  block_bits     = block_len << 3;
  terminated_len = block_bits + CCSDS_TC_STATE_BITS;

  self->length = block_bits;
  self->cw_length = terminated_len * self->params.rate;
  self->terminated_len = terminated_len;

  /* Allocate buffers */
  ALLOC(self->permutate_inv, block_bits, int);
  ALLOC(self->x1, n1[rate - 1] * terminated_len, softbit_t);
  ALLOC(self->x2, n2[rate - 1] * terminated_len, softbit_t);
  ALLOC(self->y1, terminated_len, softbit_t);
  ALLOC(self->y2, terminated_len, softbit_t);
  ALLOC(self->y0, terminated_len, softbit_t);

  memset(self->y2, 0, terminated_len * sizeof(softbit_t));
  memset(self->y0, 0, terminated_len * sizeof(softbit_t));
  
  /* Allocate permutator */
  if ((self->permutate_dir = ccsds_tc_decoder_make_permutator(block_len)) == NULL)
    goto fail;
  
  for (i = 0; i < block_bits; ++i)
    self->permutate_inv[self->permutate_dir[i]] = i;
    
  self->n1 = n1[rate - 1];

  /* Initialize encoders and finish */
  if (!ccsds_tc_bcjr_ctx_init(
    &self->dec1, 
    poly_f1, 
    poly_b, 
    n1[rate - 1], 
    block_bits)) {
    fprintf(stderr, "%s: DEC1 initialization failed\n", __FUNCTION__);
    goto fail; 
  }
  
  if (!ccsds_tc_bcjr_ctx_init(
    &self->dec2, 
    poly_f2, 
    poly_b, 
    n2[rate - 1], 
    block_bits)) {
    fprintf(stderr, "%s: DEC2 initialization failed\n", __FUNCTION__);
    goto fail; 
  }

  return true;

fail:
  ccsds_tc_decoder_finalize(self);

  return false;
}

CCSDS_TC_INLINE void
ccsds_tc_decode(ccsds_tc_decoder_t *self)
{
  unsigned int i;

  /* First permutation */
  ccsds_tc_bcjr_ctx_decode(
    &self->dec1,
    self->x1,
    self->y0,
    self->y1,
    self->params.L_c,
    self->permutate_dir);

  ccsds_tc_bcjr_ctx_decode(
    &self->dec2,
    self->x2,
    self->y1,
    self->y2,
    self->params.L_c,
    self->permutate_inv);

  /* Next permutations */
#pragma GCC unroll 2
  for (i = 1; i < self->params.iters; ++i) {
    ccsds_tc_bcjr_ctx_decode(
      &self->dec1,
      self->x1,
      self->y2,
      self->y1,
      self->params.L_c,
      self->permutate_dir);

    ccsds_tc_bcjr_ctx_decode(
      &self->dec2,
      self->x2,
      self->y1,
      self->y2,
      self->params.L_c,
      self->permutate_inv);
  }

  ccsds_tc_bcjr_invalidate_cache(&self->dec1);
  ccsds_tc_bcjr_invalidate_cache(&self->dec2);
}

void
ccsds_tc_decoder_feed_block(ccsds_tc_decoder_t *self, const softbit_t *block)
{
  (self->ser2par_sb_cb) (self, block);
  ccsds_tc_decode(self);
}

void
ccsds_tc_decoder_feed_block_float(ccsds_tc_decoder_t *self, const float *block)
{
  (self->ser2par_float_cb) (self, block);
  ccsds_tc_decode(self);
}

void
ccsds_tc_decoder_finalize(ccsds_tc_decoder_t *self)
{
  if (self->x1 != NULL)
    free(self->x1);
  
  if (self->x2 != NULL)
    free(self->x2);
  
  if (self->y0 != NULL)
    free(self->y0);

  if (self->y1 != NULL)
    free(self->y1);
  
  if (self->y2 != NULL)
    free(self->y2);

  if (self->permutate_dir != NULL)
    free(self->permutate_dir);

  if (self->permutate_inv != NULL)
    free(self->permutate_inv);

  ccsds_tc_bcjr_ctx_finalize(&self->dec1);
  ccsds_tc_bcjr_ctx_finalize(&self->dec2);
}

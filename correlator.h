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

#ifndef _CCSDS_CORRELATOR_H
#define _CCSDS_CORRELATOR_H

#include <math.h>
#include <complex.h>

#include <stdbool.h>

struct correlator {
  unsigned int seq_size;
  float *sequence;

  /* To be updated when a candidate is found */
  unsigned int where;
  float snr;
  bool reverse;
};

typedef struct correlator correlator_t;

#define correlator_INITIALIZER {0, NULL, 0, 0, 0}

bool correlator_init_from_string(correlator_t *self, const char *bits);

bool correlator_find_complex(
  correlator_t *self, 
  const float complex *x, 
  unsigned int len,
  float SNR,
  bool quadrature);

bool correlator_find(
  correlator_t *self, 
  const float *x, 
  unsigned int len, 
  float SNR);


void correlator_finalize(correlator_t *self);

static inline float
correlator_get_snr(const correlator_t *self)
{
  return self->snr;
}

static inline unsigned int
correlator_get_pos(const correlator_t *self)
{
  return self->where;
}

static inline bool
correlator_is_reverse(const correlator_t *self)
{
  return self->reverse;
}

#endif /* _CCSDS_CORRELATOR_H */

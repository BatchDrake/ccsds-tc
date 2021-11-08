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

#include <math.h>
#include <complex.h>
#include <string.h>

#include "correlator.h"
#include "defs.h"

void
correlator_finalize(correlator_t *self)
{
  if (self->sequence != NULL)
    free(self->sequence);
}

bool
correlator_init_from_string(correlator_t *self, const char *bits)
{
  unsigned int i;
  bool ok = false;

  memset(self, 0, sizeof(correlator_t));

  self->seq_size = strlen(bits);

  TRY_ALLOC(self->sequence, self->seq_size, float);

  for (i = 0; i < self->seq_size; ++i) {
    if (bits[i] == '0')
      self->sequence[i] = -1;
    else if (bits[i] == '1')
      self->sequence[i] = +1;
    else
      self->sequence[i] = 0;
  }

  ok = true;

done:
  if (!ok)
    correlator_finalize(self);

  return ok;
}

bool
correlator_find(correlator_t *self, const float *__restrict x, unsigned int len, float SNR)
{
  unsigned int i, j, p;
  float threshold = 0;
  float dotprod = 0;
  bool found = false;

  threshold = 0;
  
  if (len < self->seq_size)
    return false;

  for (i = 0; i < len; ++i)
    threshold += x[i] * x[i];

  /* 
   * This is sigma2 as seen at the correlators output in the absence
   * of correlation.
   */
  threshold *= (float) self->seq_size / (float) len;
  threshold *= SNR;

  self->snr = 0;

  for (i = 0; i < len - self->seq_size; ++i) {
    dotprod = 0;
    p = i;
    for (j = 0; j < self->seq_size; ++j)
      dotprod += self->sequence[j] * x[p++];

    if (dotprod > threshold && fabs(dotprod / threshold) > self->snr) {
      self->snr     = fabs(dotprod / threshold);
      self->reverse = dotprod < 0;
      self->where   = i;
      found = true;
    }
  }

  return found;
}

bool
correlator_find_complex(
  correlator_t *self, 
  const float complex *__restrict x, 
  unsigned int len,
  float SNR,
  bool quadrature)
{
  unsigned int i, j, p;
  float threshold = 0;
  float sigma2 = 0;
  float dotprod = 0;
  float corrpwr;
  bool found = false;

  sigma2 = 0;
  
  if (len < self->seq_size)
    return false;

  if (quadrature)
    for (i = 0; i < len; ++i)
      sigma2 += cimagf(x[i]) * cimagf(x[i]);
  else
    for (i = 0; i < len; ++i)
      sigma2 += crealf(x[i]) * crealf(x[i]);

  /* 
   * This is sigma2 as seen at the correlators output in the absence
   * of correlation.
   */
  sigma2 *= (float) self->seq_size / (float) len;
  threshold = SNR * sigma2;

  self->snr = 0;

  for (i = 0; i < len - self->seq_size; ++i) {
    dotprod = 0;
    p = i;

    if (quadrature)
      for (j = 0; j < self->seq_size; ++j)
        dotprod += self->sequence[j] * cimagf(x[p++]);
    else
      for (j = 0; j < self->seq_size; ++j)
        dotprod += self->sequence[j] * crealf(x[p++]);

    /* 
     * The idea is that the mean of the square of the correlation 
     * is given by len(seq_size) * Var(x)
     */
    corrpwr = dotprod * dotprod;

    if (corrpwr > threshold && corrpwr / sigma2 > self->snr) {
      self->snr     = corrpwr / sigma2;
      self->reverse = dotprod < 0;
      self->where   = i;
      found = true;
    }
  }

  return found;
}


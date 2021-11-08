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

#ifndef _CCSDS_GUESS_H
#define _CCSDS_GUESS_H

#include <math.h>
#include <complex.h>

#include <stdbool.h>

#define CCSDS_TC_BLOCK_SIZE_COUNT     6
#define CCSDS_TC_MINIMUM_BLOCK_LENGTH 223

struct ccsds_tc_syncword {
  const char *syncword;
  unsigned int rate_inv;
};

const struct ccsds_tc_syncword *ccsds_tc_syncword_lookup(
  unsigned int rate_inv);

/* 
 * CCSDS enables only 5 different types of block information length,
 * of which 3 is forbidden.
 */
struct ccsds_tc_report {
  const struct ccsds_tc_syncword *tc_params;
  unsigned int histogram[CCSDS_TC_BLOCK_SIZE_COUNT];
  unsigned int count_i;
  unsigned int count_q;
};

static inline unsigned int
ccsds_tc_report_count_frames(const struct ccsds_tc_report *self)
{
  return self->count_i + self->count_q;
}

static inline unsigned int
ccsds_tc_report_get_candidate_len(const struct ccsds_tc_report *self)
{
  unsigned int count, max;
  unsigned int i;

  count = 0;
  max   = 0;

  for (i = 0; i < CCSDS_TC_BLOCK_SIZE_COUNT; ++i)
    if (self->histogram[i] > count) {
      max   = i;
      count = self->histogram[i];
    }

  return max;
}

void ccsds_tc_report_init(
  struct ccsds_tc_report *self,
  const struct ccsds_tc_syncword *tcs);

bool ccsds_tc_report_find_frames(
  struct ccsds_tc_report *self,
  const float complex *x,
  unsigned int L, 
  float SNR,
  bool channel);

bool ccsds_tc_report_find_best(
  struct ccsds_tc_report *self,
  const float complex *x,
  unsigned int L,
  float SNR);

#endif /* _CCSDS_GUESS_H */


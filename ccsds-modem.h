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

#ifndef _CCSDS_MODEM_H
#define _CCSDS_MODEM_H

#include "ccsds-tc.h"
#include "correlator.h"

#include <stdint.h>

#define CCSDS_MODEM_DEFAULT_RATE      6
#define CCSDS_MODEM_DEFAULT_BLOCK_LEN CCSDS_TC_BLOCK_LENGTH_1115
#define CCSDS_MODEM_DEFAULT_ITERS     1
#define CCSDS_MODEM_DEFAULT_SYNC_SNR  80.
#define CCSDS_MODEM_DEFAULT_EBN0      0.6310
#define CCSDS_MODEM_DEFAULT_CHANNEL   true

struct ccsds_modem_params {
  unsigned int rate_inv;
  enum ccsds_tc_decoder_block_length block_len;
  unsigned int iters;
  float sync_snr;
  float EbN0;
  bool  q_channel;
};

#define ccsds_modem_params_INITIALIZER \
{                                      \
  CCSDS_MODEM_DEFAULT_RATE,            \
  CCSDS_MODEM_DEFAULT_BLOCK_LEN,       \
  CCSDS_MODEM_DEFAULT_ITERS,           \
  CCSDS_MODEM_DEFAULT_SYNC_SNR,        \
  CCSDS_MODEM_DEFAULT_EBN0,            \
  CCSDS_MODEM_DEFAULT_CHANNEL          \
}

struct ccsds_modem {
  struct ccsds_modem_params params;
  correlator_t       sync_corr;
  ccsds_tc_decoder_t decoder;
  bool               syncing;

  /* Precalculated variables */
  softbit_t       *buffer;
  float complex   *syncbuf;
  const softbit_t *y;
  unsigned int     p, q;
  unsigned int     sync_len;
  float            polarity;

  unsigned int     frame_size;
  unsigned int     block_size;
  unsigned int     enc_size;
  uint8_t         *frame_bytes;
};

typedef struct ccsds_modem ccsds_modem_t;

ccsds_modem_t *ccsds_modem_new(const struct ccsds_modem_params *params);

bool ccsds_modem_feed(
  ccsds_modem_t *self,
  const float complex *x,
  unsigned int L);

CCSDS_TC_INLINE const softbit_t *
ccsds_modem_get_frame_bits(const ccsds_modem_t *self)
{
  return self->y;
}

CCSDS_TC_INLINE unsigned int
ccsds_modem_get_block_length(const ccsds_modem_t *self)
{
  return self->block_size;
}

CCSDS_TC_INLINE float
ccsds_modem_get_error_rate(const ccsds_modem_t *self)
{
  unsigned int i, count;
  unsigned int block_size = ccsds_modem_get_block_length(self);
  bool b1, b2;

  count = 0;
  for (i = 0; i < block_size; ++i) {
    b2 = self->y[i] > 0;
    b1 = self->buffer[i * self->params.rate_inv] > 0;

    if (b1 != b2)
      ++count;
  }

  return (float) count / (float) block_size;
}

CCSDS_TC_INLINE const uint8_t *
ccsds_modem_get_frame_data(const ccsds_modem_t *self)
{
  return self->frame_bytes;
}

CCSDS_TC_INLINE unsigned int
ccsds_modem_get_frame_length(const ccsds_modem_t *self)
{
  return self->block_size >> 3;
}

void ccsds_modem_destroy(ccsds_modem_t *self);

#endif /* _CCSDS_MODEM */

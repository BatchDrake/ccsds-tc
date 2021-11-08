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


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <fcntl.h>
#include <errno.h>
#include <complex.h>
#include <netinet/in.h>

#include "ccsds-modem.h"
#include "ccsds-guess.h"

#include "lfsr.h"
#include "defs.h"

#define CCSDS_SCRAMBLER_LENGTH    255
#define CCSDS_FECF_SEED           0xffff

static bool *g_scrambler = NULL;

/* 
 * "Clarification for CCSDS CRC-16 Computation Algorithm", 
 * code by Jackson Pang, Kenneth Andres and J. Leigh Torgerson
 */

static uint16_t g_crc_table[256] = {
  0x0000, 0x1021, 0x2042, 0x3063, 0x4084, 0x50a5, 0x60c6, 0x70e7,
  0x8108, 0x9129, 0xa14a, 0xb16b, 0xc18c, 0xd1ad, 0xe1ce, 0xf1ef,
  0x1231, 0x0210, 0x3273, 0x2252, 0x52b5, 0x4294, 0x72f7, 0x62d6,
  0x9339, 0x8318, 0xb37b, 0xa35a, 0xd3bd, 0xc39c, 0xf3ff, 0xe3de,
  0x2462, 0x3443, 0x0420, 0x1401, 0x64e6, 0x74c7, 0x44a4, 0x5485,
  0xa56a, 0xb54b, 0x8528, 0x9509, 0xe5ee, 0xf5cf, 0xc5ac, 0xd58d,
  0x3653, 0x2672, 0x1611, 0x0630, 0x76d7, 0x66f6, 0x5695, 0x46b4,
  0xb75b, 0xa77a, 0x9719, 0x8738, 0xf7df, 0xe7fe, 0xd79d, 0xc7bc,
  0x48c4, 0x58e5, 0x6886, 0x78a7, 0x0840, 0x1861, 0x2802, 0x3823,
  0xc9cc, 0xd9ed, 0xe98e, 0xf9af, 0x8948, 0x9969, 0xa90a, 0xb92b,
  0x5af5, 0x4ad4, 0x7ab7, 0x6a96, 0x1a71, 0x0a50, 0x3a33, 0x2a12,
  0xdbfd, 0xcbdc, 0xfbbf, 0xeb9e, 0x9b79, 0x8b58, 0xbb3b, 0xab1a,
  0x6ca6, 0x7c87, 0x4ce4, 0x5cc5, 0x2c22, 0x3c03, 0x0c60, 0x1c41,
  0xedae, 0xfd8f, 0xcdec, 0xddcd, 0xad2a, 0xbd0b, 0x8d68, 0x9d49,
  0x7e97, 0x6eb6, 0x5ed5, 0x4ef4, 0x3e13, 0x2e32, 0x1e51, 0x0e70,
  0xff9f, 0xefbe, 0xdfdd, 0xcffc, 0xbf1b, 0xaf3a, 0x9f59, 0x8f78,
  0x9188, 0x81a9, 0xb1ca, 0xa1eb, 0xd10c, 0xc12d, 0xf14e, 0xe16f,
  0x1080, 0x00a1, 0x30c2, 0x20e3, 0x5004, 0x4025, 0x7046, 0x6067,
  0x83b9, 0x9398, 0xa3fb, 0xb3da, 0xc33d, 0xd31c, 0xe37f, 0xf35e,
  0x02b1, 0x1290, 0x22f3, 0x32d2, 0x4235, 0x5214, 0x6277, 0x7256,
  0xb5ea, 0xa5cb, 0x95a8, 0x8589, 0xf56e, 0xe54f, 0xd52c, 0xc50d,
  0x34e2, 0x24c3, 0x14a0, 0x0481, 0x7466, 0x6447, 0x5424, 0x4405,
  0xa7db, 0xb7fa, 0x8799, 0x97b8, 0xe75f, 0xf77e, 0xc71d, 0xd73c,
  0x26d3, 0x36f2, 0x0691, 0x16b0, 0x6657, 0x7676, 0x4615, 0x5634,
  0xd94c, 0xc96d, 0xf90e, 0xe92f, 0x99c8, 0x89e9, 0xb98a, 0xa9ab,
  0x5844, 0x4865, 0x7806, 0x6827, 0x18c0, 0x08e1, 0x3882, 0x28a3,
  0xcb7d, 0xdb5c, 0xeb3f, 0xfb1e, 0x8bf9, 0x9bd8, 0xabbb, 0xbb9a,
  0x4a75, 0x5a54, 0x6a37, 0x7a16, 0x0af1, 0x1ad0, 0x2ab3, 0x3a92,
  0xfd2e, 0xed0f, 0xdd6c, 0xcd4d, 0xbdaa, 0xad8b, 0x9de8, 0x8dc9,
  0x7c26, 0x6c07, 0x5c64, 0x4c45, 0x3ca2, 0x2c83, 0x1ce0, 0x0cc1,
  0xef1f, 0xff3e, 0xcf5d, 0xdf7c, 0xaf9b, 0xbfba, 0x8fd9, 0x9ff8,
  0x6e17, 0x7e36, 0x4e55, 0x5e74, 0x2e93, 0x3eb2, 0x0ed1, 0x1ef0,
};

/***********************************************************************
*
* FUNCTION:
*   ccsds_fcef_crc16
*
* INPUTS:
*   seed - (uint16_t) initial value of check bits.
*   buf - (uint8_t *) pointer to the buffer of
*   data over which you wish to generate check
*   bits.
*   len - (int) number to bytes of data in the buffer.
*
* OUTPUTS:
*
* RETURNS:
*   - (uint16_t) the checkbits.
*
* EXTERNALLY READ:
*   g_crc_table - (uint16_t)[256] the lookup table for the CCITT SDLC
*   generator polynomial (local to this module).
*
* EXTERNALLY MODIFIED:
*
* DESCRIPTION:
*   This function implements CRC generation with the CCITT SDLC error
*   polynomial (X16 + X12 + X5 + 1). You must provide it with an
*   initial seed value, a pointer to a data buffer, and the byte length
*   of the data buffer. It will return the unsigned 16-bit CRC.
*
*   You may use this function to generate a CRC over data in scattered
*   storage by making multiple calls to it. Just make sure that you
*   pass a seed of 0xFFFF on the first call. On subsequent calls, pass
*   a seed containing the return value of the previous call.
*/

CCSDS_TC_INLINE uint16_t
ccsds_fcef_crc16(uint16_t seed, const uint8_t *buf, unsigned int len)
{
  uint16_t crc = seed;
  const uint8_t *p = buf;
  
  while (len--)
    crc = g_crc_table[((crc >> 8) ^ *p++) & 0xff] ^ (crc << 8);

  return crc;
}

const bool *
ccsds_modem_get_scrambler(void)
{
  unsigned int taps[] = {0, 3, 5, 7, 8};
  unsigned int i;
  bool *result = NULL;
  lfsr_t *lfsr = NULL;

  if (g_scrambler == NULL) {
    if ((lfsr = lfsr_new(taps, 5)) == NULL)
    goto done;

    TRY_ALLOC(result, CCSDS_SCRAMBLER_LENGTH, bool);

    for (i = 0; i < CCSDS_SCRAMBLER_LENGTH; ++i)
      result[i] = lfsr_scramble(lfsr, 0);

    g_scrambler = result;  
  }
  
done:
  if (lfsr != NULL)
    lfsr_destroy(lfsr);
  
  return g_scrambler;
}

void
ccsds_modem_destroy(ccsds_modem_t *self)
{
  if (self->buffer != NULL)
    free(self->buffer);

  if (self->frame_bytes != NULL)
    free(self->frame_bytes);

  if (self->syncbuf != NULL)
    free(self->syncbuf);

  ccsds_tc_decoder_finalize(&self->decoder);
  correlator_finalize(&self->sync_corr);

  free(self);
}

ccsds_modem_t *
ccsds_modem_new(const struct ccsds_modem_params *params)
{
  ccsds_modem_t *new = NULL;
  const struct ccsds_tc_syncword *syncword;
  struct ccsds_tc_decoder_params tc_params 
    = ccsds_tc_decoder_params_INITIALIZER;

  if (ccsds_modem_get_scrambler() == NULL) {
    fprintf(stderr, "%s: failed to initialize scrambler\n", __FUNCTION__);
    return false;
  }

  if ((syncword = ccsds_tc_syncword_lookup(params->rate_inv)) == NULL) {
    fprintf(stderr, "%s: invalid rate 1/%d\n", __FUNCTION__, params->rate_inv);
    goto done;
  }

  TRY_CALLOC(new, 1, ccsds_modem_t);

  new->params = *params;

  tc_params.iters     = params->iters;
  tc_params.rate      = params->rate_inv;
  tc_params.L_c       = CCSDS_TC_EBN0_TO_LC(CCSDS_TC_F2SB(params->EbN0));
  tc_params.block_len = params->block_len;

  new->block_size = tc_params.block_len * CCSDS_TC_MINIMUM_BLOCK_LENGTH << 3;

  if (!ccsds_tc_decoder_init(&new->decoder, &tc_params))
    goto done;

  new->enc_size = ccsds_tc_decoder_get_codeword_size(&new->decoder);
  new->sync_len = strlen(syncword->syncword);
  new->syncing  = true;

  TRY_ALLOC(new->buffer,  new->enc_size, softbit_t);
  TRY_ALLOC(new->syncbuf, 2 * new->sync_len, float complex);
  TRY_ALLOC(
    new->frame_bytes, 
    tc_params.block_len * CCSDS_TC_MINIMUM_BLOCK_LENGTH,
    uint8_t);

  if (!correlator_init_from_string(&new->sync_corr, syncword->syncword)) {
    fprintf(
      stderr, 
      "%s: failed to initialize correlator for 1/%d syncword\n",
      __FUNCTION__,
      params->rate_inv);
    goto done;
  }

  return new;

done:
  if (new != NULL)
    ccsds_modem_destroy(new);

  return new;
}

CCSDS_TC_INLINE unsigned int
ccsds_modem_feed_sync(
  ccsds_modem_t *self,
  const float complex *__restrict x,
  unsigned int L)
{
  bool sync_found = false;
  unsigned int search_len = 2 * self->sync_len;
  unsigned int consumed, p, i, remainder;
  int s;

  consumed = L;
  p  = self->p;

  if (consumed > search_len - p)
    consumed = search_len - p;

  for (i = 0; i < consumed; ++i)
    self->syncbuf[p++] = x[i];

  if (p == search_len) {
    sync_found = correlator_find_complex(
      &self->sync_corr, 
      self->syncbuf, 
      search_len, 
      self->params.sync_snr,
      self->params.q_channel);
    memmove(
      self->syncbuf, 
      self->syncbuf + self->sync_len, 
      sizeof(float complex) * self->sync_len);
    p = sync_found ? 0 : self->sync_len;
  } else {
    sync_found = false;
  }

  if (sync_found) {
    /* SYNC FOUND. Populate softbits */
    self->polarity = correlator_is_reverse(&self->sync_corr) ? +1 : -1;
    /* 
      * If syncword was found in position p, we need to copy the remaining
      * search_len - (p + sync_len) softbits to the buffer. We also state
      * now that we only consumed the first p + sync_len softbits.
      */
    
    remainder = correlator_get_pos(&self->sync_corr) + self->sync_len;
    i = remainder;
    p = 0;

    /* Note that search_len is never bigger than the buffer */
    while (i < search_len) {
      s = 2 * (int) g_scrambler[p % CCSDS_SCRAMBLER_LENGTH] - 1;

      if (self->params.q_channel)
        self->buffer[p++] = s * CCSDS_TC_F2SB(
          self->polarity * cimag(self->syncbuf[i++]));
      else
        self->buffer[p++] = s * CCSDS_TC_F2SB(
          self->polarity * creal(self->syncbuf[i++]));
    }

    self->syncing = false;
  }

  self->p = p;

  return consumed;
}

CCSDS_TC_INLINE unsigned int
ccsds_modem_feed_payload(
  ccsds_modem_t *self,
  const float complex *__restrict x,
  unsigned int L)
{
  unsigned int consumed, i, p;
  int s;

  /* STATE 1: DECODING. Decide how much we should consume beforehad */
  consumed = L;
  p = self->p;

  if (consumed > self->enc_size - p)
    consumed = self->enc_size - p;

  if (self->params.q_channel) {
    for (i = 0; i < consumed; ++i, ++p) {
      s = 2 * (int) g_scrambler[p % CCSDS_SCRAMBLER_LENGTH] - 1;
      self->buffer[p] = s * CCSDS_TC_F2SB(self->polarity * cimag(x[i]));
    }
  } else {
    for (i = 0; i < consumed; ++i, ++p) {
      s = 2 * (int) g_scrambler[p % CCSDS_SCRAMBLER_LENGTH] - 1;
      self->buffer[p] = s * CCSDS_TC_F2SB(self->polarity * creal(x[i]));
    }
  }

  if (p == self->enc_size) {
    /* FULL FRAME! */
    ccsds_tc_decoder_feed_block(&self->decoder, self->buffer);
    self->y = ccsds_tc_decoder_get_output(&self->decoder);
    p = 0;
    self->syncing = true;
  }

  self->p = p;

  return consumed;
}

CCSDS_TC_INLINE bool
ccsds_modem_try_extract_frame(ccsds_modem_t *self)
{
  unsigned int i;
  uint16_t crc;
  unsigned int frame_len = self->block_size >> 3;

  for (i = 0; i < self->block_size; ++i)
    if ((i & 7) == 0)
      self->frame_bytes[i >> 3] = (self->y[i] > 0) << 7;
    else
      self->frame_bytes[i >> 3] |= (self->y[i] > 0) << (7 - (i & 7));

  /* CCSDS frames are big endian */
  crc = ntohs(*(const uint16_t *) &self->frame_bytes[frame_len - 2]);

  /* Check FECF */
  return ccsds_fcef_crc16(
    CCSDS_FECF_SEED, 
    self->frame_bytes, 
    frame_len - 2) == crc;
}

bool 
ccsds_modem_feed(
  ccsds_modem_t *self,
  const float complex *__restrict x,
  unsigned int L)
{
  bool have_frame = false;
  unsigned int consumed;

  while (L > 0) {
    if (self->syncing) {
      consumed = ccsds_modem_feed_sync(self, x, L);
    } else {
      consumed = ccsds_modem_feed_payload(self, x, L);
      have_frame = have_frame || self->syncing;
    }

    x += consumed;
    L -= consumed;
  }

  /* Something that looks like a frame? */
  if (have_frame)
    have_frame = ccsds_modem_try_extract_frame(self);

  return have_frame;
}

bool 
ccsds_modem_process_area(
  ccsds_modem_t *self,
  const float complex *__restrict x,
  unsigned int L)
{
  float polarity = 1, f;
  const bool *scrambler = NULL;
  unsigned int p = 0, i, q;
  unsigned int block_size;
  unsigned int enc_size;
  unsigned int sync_step;
  unsigned int frame_size;
  bool syncing = true;
  bool ok = false;
 
  enc_size   = self->enc_size;
  sync_step  = self->sync_len;
  frame_size = sync_step + enc_size + 4 * self->params.rate_inv;

  if ((scrambler = ccsds_modem_get_scrambler()) == NULL)
    goto done;

  /* Decoding loop */
  while (p < L - frame_size) {
    if (syncing) {
      if (correlator_find_complex(
        &self->sync_corr, 
        x + p, 
        2 * sync_step, 
        self->params.sync_snr, 
        self->params.q_channel)) {
        polarity = correlator_is_reverse(&self->sync_corr) ? -1 : +1;
        p += correlator_get_pos(&self->sync_corr);
        syncing  = false;
      }

      p += sync_step;
    } else {
      q = 0;
      if (self->params.q_channel) {
        for (i = 0; i < enc_size; ++i) {
          f = 2 * (int) scrambler[q++ % CCSDS_SCRAMBLER_LENGTH] - 1;
          self->buffer[i] = CCSDS_TC_F2SB(-polarity * cimag(x[p + i]) * f);
        }
      } else {
        for (i = 0; i < enc_size; ++i) {
          f = 2 * (int) scrambler[q++ % CCSDS_SCRAMBLER_LENGTH] - 1;
          self->buffer[i] = CCSDS_TC_F2SB(-polarity * creal(x[p + i]) * f);
        }
      }

      ccsds_tc_decoder_feed_block(&self->decoder, self->buffer);
      self->y = ccsds_tc_decoder_get_output(&self->decoder);

      /* TODO: Trigger action */    
      p += enc_size;
      syncing = true;
    }
  }

  ok = true;

done:
  return ok;
}


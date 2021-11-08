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

#include "ccsds-guess.h"
#include "correlator.h"

#include <stdio.h>
#include <string.h>

#define CCSDS_TC_SYNC_RATE_1_2 "0000001101000111011101101100011100100111001010001001010110110000"
#define CCSDS_TC_SYNC_RATE_1_3 "001001011101010111000000110011101000100110010000111101101100100101000110000110111111011110011100"
#define CCSDS_TC_SYNC_RATE_1_4 "00000011010001110111011011000111001001110010100010010101101100001111110010111000100010010011100011011000110101110110101001001111"
#define CCSDS_TC_SYNC_RATE_1_6 "001001011101010111000000110011101000100110010000111101101100100101000110000110111111011110011100110110100010101000111111001100010111011001101111000010010011011010111001111001000000100001100011"

#define CCSDS_TC_SYNCWORD_COUNT 4

static const struct ccsds_tc_syncword g_syncwords[] = {
  {CCSDS_TC_SYNC_RATE_1_2, 2},
  {CCSDS_TC_SYNC_RATE_1_3, 3},
  {CCSDS_TC_SYNC_RATE_1_4, 4},
  {CCSDS_TC_SYNC_RATE_1_6, 6}
};

void
ccsds_tc_report_init(struct ccsds_tc_report *self, const struct ccsds_tc_syncword *tcs)
{
  memset(self, 0, sizeof(struct ccsds_tc_report));
  self->tc_params = tcs;
}

const struct ccsds_tc_syncword *
ccsds_tc_syncword_lookup(unsigned int rate_inv)
{
  unsigned int i;

  for (i = 0; i < 4; ++i)
    if (g_syncwords[i].rate_inv == rate_inv)
      return &g_syncwords[i];

  return NULL;
}

static inline unsigned int
ccsds_tc_report_get_block_length(const struct ccsds_tc_report *self, unsigned int L)
{
  unsigned int payload;
  unsigned int len;
  unsigned int rateinv = self->tc_params->rate_inv;

  if (L % rateinv != 0)
    return 0;

  L /= rateinv;

  payload = L - 32 - 4;

  if (payload % (CCSDS_TC_MINIMUM_BLOCK_LENGTH << 3) == 0) {
    len = payload / (CCSDS_TC_MINIMUM_BLOCK_LENGTH << 3);
    if (len < 1 || len > 5)
      return 0;

    return len;
  }

  return 0;
}

bool
ccsds_tc_report_find_frames(
  struct ccsds_tc_report *self,
  const float complex *__restrict x,
  unsigned int L, 
  float SNR,
  bool channel)
{
  correlator_t correlator = correlator_INITIALIZER;
  unsigned int i;
  unsigned int curr, last;
  unsigned int step;
  unsigned int block_info_len;

  bool ok = false;

  step = strlen(self->tc_params->syncword);

  if (L < 2 * step)
    return true;

  if (!correlator_init_from_string(&correlator, self->tc_params->syncword)) {
    fprintf(
      stderr, 
      "%s: failed to initialize correlator for 1/%d syncword\n",
      __FUNCTION__,
      self->tc_params->rate_inv);
    goto done;
  }

  last = L;
  
  for (i = 0; i < L - 2 * step; i += step) {
    if (correlator_find_complex(&correlator, x + i, 2 * step, SNR, channel)) {
      curr = correlator_get_pos(&correlator) + i;

      if (
        last < curr 
        && (block_info_len = ccsds_tc_report_get_block_length(self, curr - last)) > 0) {
        printf(
        "%c channel: Potential 1/%d CCSDS frame @ pos = %7d (%7d sice last) (SNR: %.1f dB)\r",
        channel ? 'Q' : 'I',
        self->tc_params->rate_inv,
        correlator_get_pos(&correlator) + i,
        curr - last,
        10 * log10(correlator_get_snr(&correlator)));

        ++self->histogram[block_info_len];
        
        if (channel)
          ++self->count_q;
        else
          ++self->count_i;

        fflush(stdout);
      }

      last = curr;
    }
  }

  if ((self->count_i + self->count_q) > 0)
    printf("\033[2K");

  ok = true;

done:
  correlator_finalize(&correlator);

  return ok;
}

bool
ccsds_tc_report_find_best(
  struct ccsds_tc_report *self,
  const float complex *__restrict x,
  unsigned int L,
  float SNR)
{
  struct ccsds_tc_report report, best_report;
  unsigned int i;
  bool ok = false;

  memset(&best_report, 0, sizeof(struct ccsds_tc_report));

  for (i = 0; i < CCSDS_TC_SYNCWORD_COUNT; ++i) {
    printf("Looking for syncwords for r = 1/%d\n", g_syncwords[i].rate_inv);
    ccsds_tc_report_init(&report, &g_syncwords[i]);

    ccsds_tc_report_find_frames(&report, x, L, SNR, false);
    ccsds_tc_report_find_frames(&report, x, L, SNR, true);

    if (ccsds_tc_report_count_frames(&report) 
      > ccsds_tc_report_count_frames(&best_report))
      best_report = report;
  }

  if (ccsds_tc_report_count_frames(&best_report) > 0) {
    printf("Candidate CCSDS turbocode found\n");
    printf("    Code rate:    1/%d\n", best_report.tc_params->rate_inv);
    if (best_report.count_i > best_report.count_q)
      printf("    Channel I:    %d frames\n", best_report.count_i);
    else
      printf("    Channel Q:    %d frames\n", best_report.count_q);

    printf(
      "    Frame length: %d bits\n", 
      ccsds_tc_report_get_candidate_len(&best_report) 
      * CCSDS_TC_MINIMUM_BLOCK_LENGTH << 3);
  } else {
    printf("No valid CCSDS frames found\n");
  }

  *self = best_report;

  ok = true;

done:
  return ok;
}

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
#include <stdlib.h>
#include <sys/mman.h>

#include "ccsds-tc.h"

void
parse_packet(const softbit_t *u, const float *r)
{
  ccsds_tc_decoder_t decoder;
  const softbit_t *y;
  struct ccsds_tc_decoder_params params = ccsds_tc_decoder_params_INITIALIZER;
  unsigned int i;
  unsigned int length;
  bool b1, b2, b3;

  params.block_len = CCSDS_TC_BLOCK_LENGTH_223;
  params.iters     = 2;
  params.L_c       = CCSDS_TC_EBN0_TO_LC(CCSDS_TC_F2SB(powf(10, .1 * -2)));
  params.rate      = 6;
  
  if (!ccsds_tc_decoder_init(&decoder, &params)) {
    fprintf(stderr, "%s: failed to create decoder\n", __FUNCTION__);
    goto done;
  }

  ccsds_tc_decoder_feed_block_float(&decoder, r);

  length = ccsds_tc_decoder_get_block_size(&decoder);
  y = ccsds_tc_get_output(&decoder);

  for (i = 0; i < length; ++i) {
    b1 = u[i] > 0;
    b2 = r[i * 6] > 0;
    b3 = y[i] > 0;

    if (b1 != b2 || b2 != b3) {
      printf("BIT %5d / %5d CORRUPTED: ", i + 1, length);
      if (b1 == b3) {
        printf("\033[1;32mCORRECTED\033[0m\n");
      } else {
        printf("\033[1;31mFAILED (L(u) = %g)\033[0m\n", CCSDS_TC_SB2F(y[i]));
      }
    }
  }

done:
  ccsds_tc_decoder_finalize(&decoder);
}

int
main(int argc, char **argv)
{
  size_t size;
  void *map, *map2;
  FILE *fp = fopen("turbosignal.raw", "rb");
  FILE *fp2 = fopen("rightsignal.raw", "rb");

  if (fp == NULL) {
    perror("turbosignal.raw");
    return 1;
  }
  
  if (fp2 == NULL) {
    perror("rightsignal.raw");
    return 1;
  }

  fseek(fp, 0, SEEK_END);
  size = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  size /= sizeof(float);

  map  = mmap(NULL, size * sizeof(float), PROT_READ, MAP_PRIVATE, fileno(fp), 0);
  if (map == (caddr_t) - 1) {
    perror("mmap turbosignal");
    exit(EXIT_FAILURE);
  }

  map2 = mmap(NULL, size * sizeof(float) / 6, PROT_READ, MAP_PRIVATE, fileno(fp2), 0);
  if (map == (caddr_t) - 1) {
    perror("mmap rightsignal");
    exit(EXIT_FAILURE);
  }
  fclose(fp);

  parse_packet(map2, map);

  return 0;
}
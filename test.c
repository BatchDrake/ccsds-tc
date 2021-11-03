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
#include <string.h>
#include <sys/mman.h>
#include <sys/time.h>

#include "ccsds-tc.h"

#define CCSDS_TC_TEST_ITERS 3000

enum ccsds_tc_test_type {
  CCSDS_TC_TEST_TYPE_SINGLE_RUN,
  CCSDS_TC_TEST_TYPE_DECODER_RATE,
  CCSDS_TC_TEST_TYPE_ITERATION_RATE
};

void
run_benchmark(const softbit_t *u, const float *r, enum ccsds_tc_test_type type)
{
  ccsds_tc_decoder_t decoder;
  const softbit_t *y;
  float rate;
  struct ccsds_tc_decoder_params params = ccsds_tc_decoder_params_INITIALIZER;
  struct timeval tv, otv, diff;
  unsigned int i;
  unsigned int length;
  bool b1, b2, b3;

  params.block_len = CCSDS_TC_BLOCK_LENGTH_223;
  params.iters     = type == CCSDS_TC_TEST_TYPE_ITERATION_RATE 
    ? CCSDS_TC_TEST_ITERS
    : 1;

  params.L_c       = CCSDS_TC_EBN0_TO_LC(CCSDS_TC_F2SB(powf(10, .1 * -2)));
  params.rate      = 6;
  
  if (!ccsds_tc_decoder_init(&decoder, &params)) {
    fprintf(stderr, "%s: failed to create decoder\n", __FUNCTION__);
    goto done;
  }

  length = ccsds_tc_decoder_get_block_size(&decoder);

  switch (type) {
    case CCSDS_TC_TEST_TYPE_SINGLE_RUN:
      printf("Running one-shot decoding with %d iters\n", params.iters);
      ccsds_tc_decoder_feed_block_float(&decoder, r);
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

      break;

    case CCSDS_TC_TEST_TYPE_DECODER_RATE:
      printf(
        "Decoding the same %d-bit block %d times, press Ctrl+C to abort\n",
        ccsds_tc_decoder_get_block_size(&decoder),
        CCSDS_TC_TEST_ITERS);
#ifdef CCSDS_TC_INT_ARITHMETICS
      printf("Note: integer arithmetics enabled. Expect overhead due to float-to-softbit\n");
      printf("conversions prior to each decoding.\n");
#endif /* CCSDS_TC_INT_ARITHMETICS */

      for (;;) {
        gettimeofday(&otv, NULL);
        for (i = 0; i < CCSDS_TC_TEST_ITERS; ++i)
          ccsds_tc_decoder_feed_block_float(&decoder, r);
        gettimeofday(&tv, NULL);

        gettimeofday(&tv, NULL);
        timersub(&tv, &otv, &diff);

        rate = (CCSDS_TC_TEST_ITERS * length) / (diff.tv_sec + 1e-6 * diff.tv_usec);

        printf("Output rate: ");
        if (rate < 1e3)
          printf("%.3f bps\n", rate);
        else if (rate < 1e6)
          printf("%.3f kbps\n", rate * 1e-3);
        else
          printf("%.3f Mbps\n", rate * 1e-6);
      }
      break;

    case CCSDS_TC_TEST_TYPE_ITERATION_RATE:
      printf(
        "Performing %d iterations over the same %d-bit block, press Ctrl+C to abort\n",
        CCSDS_TC_TEST_ITERS,
        ccsds_tc_decoder_get_block_size(&decoder));

      for (;;) {
        gettimeofday(&otv, NULL);
        ccsds_tc_decoder_feed_block_float(&decoder, r);
        gettimeofday(&tv, NULL);

        gettimeofday(&tv, NULL);
        timersub(&tv, &otv, &diff);

        rate = (CCSDS_TC_TEST_ITERS * length) / (diff.tv_sec + 1e-6 * diff.tv_usec);

        printf("Iteration rate: ");
        if (rate < 1e3)
          printf("%.3f bps\n", rate);
        else if (rate < 1e6)
          printf("%.3f kbps\n", rate * 1e-3);
        else
          printf("%.3f Mbps\n", rate * 1e-6);
      }
      break;
  }
  
done:
  ccsds_tc_decoder_finalize(&decoder);
}

int
main(int argc, char **argv)
{
  size_t size;
  enum ccsds_tc_test_type type = CCSDS_TC_TEST_TYPE_SINGLE_RUN;
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

  if (argc < 2) {
    fprintf(stderr, "%s: no benchmark specified, defaulting to single\n", argv[0]);
  } else {
    if (strcmp(argv[1], "single") == 0)
      type = CCSDS_TC_TEST_TYPE_SINGLE_RUN;
    else if (strcmp(argv[1], "rate") == 0)
      type = CCSDS_TC_TEST_TYPE_DECODER_RATE;
    else if (strcmp(argv[1], "iter") == 0)
      type = CCSDS_TC_TEST_TYPE_ITERATION_RATE;
    else {
      fprintf(stderr, "%s: invalid benchmark type `%s'\n", argv[0], argv[1]);
      exit(EXIT_FAILURE);
    }
  }

  run_benchmark(map2, map, type);

  return 0;
}

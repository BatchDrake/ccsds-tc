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

#ifndef _LFSR_H
#define _LFSR_H

#include <stdint.h>

#define LFSR_MAX_TAPS 63

struct lfsr {
  uint64_t mask;
  uint64_t reg;
  uint64_t len;
  uint64_t cycle_len; /* Assuming it's primitive */
};

typedef struct lfsr lfsr_t;

static inline uint64_t
lfsr_get_cycle_len(const lfsr_t *self)
{
  return self->cycle_len;
}

lfsr_t *lfsr_new(const unsigned int *taps, unsigned int tap_len);
uint8_t lfsr_scramble(lfsr_t *self, uint8_t input);
uint8_t lfsr_descramble(lfsr_t *self, uint8_t input);
void lfsr_reset(lfsr_t *self);
char *lfsr_get_poly(const lfsr_t *self);
void lfsr_destroy(lfsr_t *);

#endif /* _LFSR_H */

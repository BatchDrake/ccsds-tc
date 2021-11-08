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
#include <stdarg.h>

#include "lfsr.h"

#define TRY(expr) if (!(expr)) goto fail

#define STRBUILD_BSIZ           16

static inline unsigned int
popcount64(uint64_t b)
{
  b = (b & 0x5555555555555555ull) + (b >> 1 & 0x5555555555555555ull);
  b = (b & 0x3333333333333333ull) + (b >> 2 & 0x3333333333333333ull);
  b = b + (b >> 4) & 0x0F0F0F0F0F0F0F0Full;
  b = b + (b >> 8);
  b = b + (b >> 16);
  b = b + (b >> 32) & 0x0000007F;

  return (unsigned int) b;
}

int
is_asciiz(const char *buf, int lbound, int ubound)
{
  register int i;

  for (i = lbound; i < ubound; i++)
    if (!buf[i])
      return i + 1;
  return 0;
}

char *
vstrbuild(const char *fmt, va_list ap)
{
  char *out;
  int size, zeroindex;
  int last;
  va_list copy;
  
  last = 0;
  
  if (fmt != NULL) {
    if (!*fmt) {
      out = malloc(1);
      out[0] = '\0';
      return out;
    }
    
    va_copy(copy, ap);
    size = vsnprintf(NULL, 0, fmt, copy) + 1;
    va_end(copy);
    
    if ((out = malloc(size)) == NULL)
      return NULL;
    
    va_copy(copy, ap);
    vsnprintf(out, size, fmt, copy);
    va_end(copy);
    
    for (;;) {
      if ((zeroindex = is_asciiz(out, last, size)) != 0)
        break;

      last = size;
      size += STRBUILD_BSIZ;
      
      out = realloc(out, size); /* Reasignamos */
      
      va_copy (copy, ap);
      vsnprintf(out, size, fmt, copy);
      va_end (copy);
    }
  }
  else
    out = NULL;
  
  return out;
}


/* Construye una cadena mediante el formato printf y devuelve un
   puntero a la cadena resultado. DEBES liberar tu mismo la salida. */

/* FIXME: Buscar alguna alternativa mas portable */
char*
strbuild (const char *fmt, ...)
{
  char *out;
  va_list ap;

  va_start (ap, fmt);
  out =  vstrbuild (fmt, ap);
  va_end (ap);

  return out;
}

void
lfsr_destroy(lfsr_t *self)
{
  free(self);
}

char *
lfsr_get_poly(const lfsr_t *self)
{
  unsigned int i;
  char *prev = NULL;
  char *poly = NULL;

  for (i = 63; i >= 1; --i)
    if ((self->mask & (1ull << i)) != 0) {
      TRY(poly = strbuild("%sx^%d + ", prev == NULL ? "" : prev, i));
      if (prev != NULL)
        free(prev);
      prev = poly;
    }

  TRY(poly = strbuild("%s1", prev == NULL ? "" : prev));

  if (prev != NULL)
    free(prev);

  return poly;

fail:
  if (prev != NULL)
    free(prev);

  return NULL;
}

void
lfsr_reset(lfsr_t *self)
{
  self->reg = (1ull << (2 * self->len)) - 1;
}

lfsr_t *
lfsr_new(const unsigned int *taps, unsigned int tap_len)
{
  lfsr_t *self = NULL;
  unsigned int i;

  if ((self = calloc(1, sizeof(lfsr_t))) == NULL)
    goto fail;

  for (i = 0; i < tap_len; ++i) {
    if (taps[i] >= LFSR_MAX_TAPS) {
      fprintf(stderr, "Invalid tap %d\n", taps[i]);
      goto fail;
    }

    self->mask |= 1ull << taps[i];
    if (self->len < taps[i])
      self->len = taps[i];
  }

  lfsr_reset(self);

  self->cycle_len = (1ull << self->len) - 1;

  --self->len;

  return self;

fail:
  if (self != NULL)
    lfsr_destroy(self);

  return NULL;
}

static inline uint8_t
lfsr_core(lfsr_t *self, uint8_t input, unsigned int direction)
{
  uint8_t x = input & 1;
  uint8_t y = (popcount64(self->reg & self->mask) & 1) ^ x;
  uint8_t newbit = direction ? x : y;
  uint8_t output = direction ? y : self->reg & 1;

  self->reg = (self->reg >> 1) | (newbit << self->len);

  return y;
}

uint8_t
lfsr_scramble(lfsr_t *self, uint8_t input)
{
  return lfsr_core(self, input, 0);
}

uint8_t
lfsr_descramble(lfsr_t *self, uint8_t input)
{
  return lfsr_core(self, input, 1);
}


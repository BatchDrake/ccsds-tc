#ifndef _DEFS_H
#define _DEFS_H

#include <stdio.h>
#include <stdlib.h>

#define TRY_ALLOC(where, N, what)                         \
  if ((where = malloc(N * sizeof(what))) == NULL) {       \
    fprintf(                                              \
      stderr,                                             \
      "%s: failed to allocate %d elements of type %s\n",  \
      __FILE__,                                           \
      N,                                                  \
      #what);                                             \
    goto done;                                            \
  }

#define TRY_CALLOC(where, N, what)                        \
  if ((where = calloc(N, sizeof(what))) == NULL) {        \
    fprintf(                                              \
      stderr,                                             \
      "%s: failed to allocate %d elements of type %s\n",  \
      __FILE__,                                           \
      N,                                                  \
      #what);                                             \
    goto done;                                            \
  }

#ifndef __GNUC__
#  ifdef __restricted
#    undef __restricted
#  endif /* __restricted */
#  define __restricted
#endif /* __GNUC__ */

#endif /* _DEFS_H */

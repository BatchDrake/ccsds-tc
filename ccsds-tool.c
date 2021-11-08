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
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>
#include <complex.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <signal.h>
#include <netdb.h>

#include "ccsds-modem.h"
#include "ccsds-guess.h"
#include "defs.h"

#define CCSDS_TOOL_VERSION                "0.1"
#define CCSDS_TOOL_RATE_REPORT_INTERVAL_S 1
#define CCSDS_TOOL_READ_BUFFER_SIZE       (4096 / sizeof(float complex))

struct ccsds_tool_options {
  const char *guess;
  const char *file;
  const char *ofile;
  unsigned int port;
  unsigned int rate_inv;
  bool q_channel;
  float snr;
  unsigned int iters;
  unsigned int block_size;
};

#define ccsds_tool_options_INITIALIZER  \
{                                       \
  NULL,  /* guess */                    \
  NULL,  /* file */                     \
  NULL,  /* ofile */                    \
  0,     /* port */                     \
  6,     /* rate_inv */                 \
  true,  /* q_channel */                \
  19.03, /* snr */                      \
  1,     /* iters */                    \
  8920,  /* block_size */               \
}

static int
open_listening_socket(const char *argv0, uint16_t port)
{
  struct sockaddr_in server_addr;
  struct sockaddr_in client_addr;
  unsigned int structsiz;
  int sfd, cfd = -1;
  int flg;
  
  flg = 1;
  
  structsiz = sizeof(struct sockaddr_in);
  
  server_addr.sin_family      = AF_INET;
  server_addr.sin_addr.s_addr = INADDR_ANY;
  server_addr.sin_port        = htons(port);
  
  if ((sfd = socket(PF_INET, SOCK_STREAM, 0)) == -1) {
    fprintf(
      stderr, 
      "%s: cannot create socket: %s\n", 
      argv0, 
      strerror(errno));
    goto done;
  }
  
  if (setsockopt(
    sfd, 
    SOL_SOCKET, 
    SO_REUSEADDR, 
    &flg, 
    sizeof(int))) {
    fprintf(
      stderr, 
      "%s: cannot set SO_REUSEADDR: %s\n", 
      argv0, 
      strerror(errno));
    goto done;
  }

  if (bind(
    sfd, 
    (struct sockaddr*) &server_addr, 
    sizeof (struct sockaddr_in))) {
    fprintf(
      stderr, 
      "%s: cannot bind socket address: %s\n", 
      argv0, 
      strerror(errno));
    goto done;
  }

  if (listen(sfd, 5) == -1) {
    fprintf(
      stderr, 
      "%s: cannot listen on socket: %s\n", 
      argv0, 
      strerror(errno));
    goto done;
  }
  
  signal(SIGPIPE, SIG_IGN);
    
  cfd = accept(
    sfd, 
    (struct sockaddr *) &client_addr,
    &structsiz);
  
  if (cfd == -1) {
    fprintf(
      stderr, 
      "%s: cannot accept connection: %s\n",
      argv0,
      strerror(errno));
    goto done;
  }

done:
  if (sfd != -1)
    close(sfd);

  return cfd;
}

CCSDS_TC_INLINE void
report_rate(unsigned int count, const struct timeval *diff)
{
  float rate = (count << 3) / (diff->tv_sec + diff->tv_usec * 1e-6);

  if (rate < 1e3)
    fprintf(stderr, "Decode rate: %.2f bps\r", rate);
  else if (rate < 1e6)
    fprintf(stderr, "Decode rate: %.2f kbps\r", rate * 1e-3);
  else
    fprintf(stderr, "Decode rate: %.2f Mbps\r", rate * 1e-6);
}

bool
run_guess(const char *argv0, const struct ccsds_tool_options *opt)
{
  struct ccsds_tc_report report;
  struct stat sbuf;
  unsigned int L;
  int fd = -1;
  float complex *samples = (float complex *) -1;
  bool ok = false;

  if (stat(opt->guess, &sbuf) == -1) {
    fprintf(
      stderr, 
      "%s: cannot stat `%s': %s\n", 
      argv0, 
      opt->guess, 
      strerror(errno));
    exit(EXIT_FAILURE);  
  }

  if ((fd = open(opt->guess, O_RDONLY)) == -1) {
    fprintf(
      stderr, 
      "%s: cannot open `%s' for reading: %s\n", 
      argv0, 
      opt->guess, 
      strerror(errno));
    exit(EXIT_FAILURE);
  }

  if ((samples = mmap(
    NULL, 
    sbuf.st_size, 
    PROT_READ, 
    MAP_PRIVATE, 
    fd, 
    0)) == (const float complex *) -1) {
    fprintf(
      stderr, 
      "%s: cannot mmap `%s': %s\n", 
      argv0, 
      opt->guess, 
      strerror(errno));
    exit(EXIT_FAILURE);
  }

  L = sbuf.st_size / sizeof(float complex);

  if (!ccsds_tc_report_find_best(
    &report, 
    samples, 
    L, 
    pow(10, .1 * opt->snr)))
    goto done;
    
  ok = true;

done:
  if (fd != -1)
    close(fd);

  if (samples != (const float complex *) -1)
    munmap(samples, sbuf.st_size);

  return ok;
}

static bool
run_decoder(const char *argv0, const struct ccsds_tool_options *opt)
{
  FILE *ifp = stdin;
  FILE *ofp = stdout;
  float complex *buffer = NULL;
  int sfd = -1;
  size_t got;
  unsigned int count = 0, total = 0;
  ccsds_modem_t *modem = NULL;
  struct ccsds_modem_params params = ccsds_modem_params_INITIALIZER;
  float rate;
  struct timeval tv, otv, diff;  
  bool ok = false;

  params.block_len = opt->block_size / (CCSDS_TC_MINIMUM_BLOCK_LENGTH << 3);
  params.iters     = opt->iters;
  params.q_channel = opt->q_channel;
  params.sync_snr  = powf(10, .1 * opt->snr);

  fprintf(stderr, "CCSDS tool v" CCSDS_TOOL_VERSION " for the Amateur DSN by EA1IYR\n");
  fprintf(stderr, "(c) 2021 Gonzalo J. Carracedo - https://actinid.org\n");
  fprintf(stderr, "  Code rate:       1/%d\n", opt->rate_inv);
  fprintf(stderr, "  Block length:    %d bits\n", opt->block_size);
  fprintf(stderr, "  Turbocode iters: %d\n", params.iters);
  fprintf(stderr, "  Channel:         %c\n", opt->q_channel ? 'Q' : 'I');
  fprintf(stderr, "  Sync SNR:        %g dB\n", opt->snr);

  TRY_ALLOC(buffer, CCSDS_TOOL_READ_BUFFER_SIZE, float complex);

  if ((modem = ccsds_modem_new(&params)) == NULL)
    goto done;

  /* Open input file */
  if (opt->port != 0) {
    fprintf(stderr, "  Waiting for connections to port %u (TCP)\n", opt->port);
    if ((sfd = open_listening_socket(argv0, opt->port)) == -1)
      goto done;
    if ((ifp = fdopen(sfd, "rb")) == NULL) {
      fprintf(
        stderr, 
        "%s: failed to open socket as file: %s\n", 
        argv0, 
        strerror(errno));
      goto done;
    }
  } else if (opt->file != NULL) {
    if ((ifp = fopen(opt->file, "rb")) == NULL) {
      fprintf(
        stderr, 
        "%s: failed to open `%s' for reading: %s\n", 
        argv0, 
        strerror(errno));
      goto done;
    }
    fprintf(stderr, "  Input file: %s\n", opt->file);
  } else {
    fprintf(stderr, "  Input file: stdin\n");
  }

  /* Open output file */
  if (opt->ofile) {
    if ((ofp = fopen(opt->ofile, "wb")) == NULL) {
      fprintf(
        stderr, 
        "%s: failed to open `%s' for writing: %s\n", 
        argv0, 
        strerror(errno));
      goto done;
    }
    fprintf(stderr, "  Output file: %s\n", opt->ofile);
  } else {
    fprintf(stderr, "  Output file: stdout\n");
  }

  /* Reader loop */
  gettimeofday(&otv, NULL);
  while ((got = 
    fread(
      buffer, 
      sizeof(float complex), 
      CCSDS_TOOL_READ_BUFFER_SIZE, 
      ifp)) > 0) {
    if (ccsds_modem_feed(modem, buffer, got)) {
      /* Frame found with valid CRC! */
      fwrite(
        ccsds_modem_get_frame_data(modem),
        ccsds_modem_get_frame_length(modem),
        1,
        ofp);
      count += ccsds_modem_get_frame_length(modem);
      total += ccsds_modem_get_frame_length(modem);
      gettimeofday(&tv, NULL);
      timersub(&tv, &otv, &diff);

      if (diff.tv_sec >= CCSDS_TOOL_RATE_REPORT_INTERVAL_S) {
        report_rate(count, &diff);
        otv = tv;
        count = 0;
      }
    }
  }
  
  gettimeofday(&tv, NULL);
  timersub(&tv, &otv, &diff);
  report_rate(count, &diff);
  fputc('\n', stderr);
  fprintf(stderr, "%s: %d bytes decoded\n", argv0, total);

  ok = true;

done:
  if (modem != NULL)
    ccsds_modem_destroy(modem);

  if (sfd != -1)
    close(sfd);

  if (ofp != NULL && ofp != stdout)
    fclose(ofp);

  if (ifp != NULL && ifp != stdin)
    fclose(ifp);

  if (buffer != NULL)
    free(buffer);

  return ok;
}

static bool
run(const char *argv0, const struct ccsds_tool_options *opt)
{
  bool ok;

  if (opt->guess != NULL)
    ok = run_guess(argv0, opt);
  else
    ok = run_decoder(argv0, opt);

  return ok;
}

void
help(const char *argv0)
{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "    %s [options]\n\n", argv0);
  fprintf(stderr, "Attempts to guess and decode CSSDS turbo-coded frames from demodulated\n");
  fprintf(stderr, "32-bit complex float samples.\n\n");

  fprintf(stderr, "OPTIONS:\n");
  fprintf(stderr, "    -g, --guess file     Attempt to guess turbocode parameters from input file\n");
  fprintf(stderr, "    -r, --rate 1/k       Set the code rate to 1/k (default is 1/6)\n");
  fprintf(stderr, "    -s, --block-size LEN Set the block size to LEN bits (default is 8920)\n");
  fprintf(stderr, "    -S, --sync-snr LEVEL Set the syncword correlation SNR to LEVEL dB\n");
  fprintf(stderr, "                         (default is 19.03 dB)\n");
  fprintf(stderr, "    -c, --component C    Sets the complex channel to use (I or Q, default is\n");
  fprintf(stderr, "                         Q: quadrature channel)\n");
  fprintf(stderr, "    -i, --iters NUM      Set the number of iterations of the BCJR turbo code\n");
  fprintf(stderr, "                         decoder (default is 1)\n");
  fprintf(stderr, "    -f, --file PATH      Read complex samples from a file instead of the\n");
  fprintf(stderr, "                         standard input\n");
  fprintf(stderr, "    -o, --output PATH    Write decoded frames to PATH instead of the\n");
  fprintf(stderr, "                         standard output\n");
  fprintf(stderr, "    -l, --listen PORT    Read complex samples from a listening socket at port\n");
  fprintf(stderr, "                         PORT instead of the standard input\n\n");
  fprintf(stderr, "    -h, --help           This help\n\n");

  fprintf(stderr, "Copyright (C) 2021 Gonzalo J. Carracedo (EA1IYR)\n");
  fprintf(stderr, "Hints on BCJR optimizations by r00t and Daniel Estévez (EA4GPZ)\n\n");
  fprintf(stderr, "This is free software; see the source for copying conditions.  There is NO\n");
  fprintf(stderr, "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
}

static struct option g_options[] = {
    {"rate",       1, 0, 'r'},
    {"block-size", 1, 0, 's'},
    {"sync-snr",   1, 0, 'S'},
    {"component",  1, 0, 'c'},
    {"iters",      1, 0, 'i'},
    {"file",       1, 0, 'f'},
    {"output",     1, 0, 'o'},
    {"listen",     1, 0, 'l'},
    {"guess",      1, 0, 'g'},
    {"help",       0, 0, 'h'},
    {0,            0, 0, 0}
};

int
main(int argc, char **argv)
{
  int c;
  int digit_optind = 0;
  int option_index;
  struct ccsds_tool_options options = ccsds_tool_options_INITIALIZER;
  
  while ((c = getopt_long(
      argc, 
      argv, 
      "r:s:S:c:i:g:f:o:l:h",
      g_options,
      &option_index)) != -1) {

    switch (c) {
      case 'r':
        if (sscanf(optarg, "1/%u", &options.rate_inv) != 1) {
          fprintf(stderr, "%s: invalid rate `%s'\n", argv[0], optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 's':
        if (sscanf(optarg, "%u", &options.block_size) != 1) {
          fprintf(stderr, "%s: invalid block size `%s'\n", argv[0], optarg);
          exit(EXIT_FAILURE);
        }

        if ((options.block_size % (CCSDS_TC_MINIMUM_BLOCK_LENGTH << 3)) != 0) {
          fprintf(
            stderr, 
            "%s: invalid block size %d: it must be a multiple of %d\n", 
            argv[0],
            options.block_size,
            CCSDS_TC_MINIMUM_BLOCK_LENGTH << 3);
          exit(EXIT_FAILURE);
        }
        break;

      case 'S':
        if (sscanf(optarg, "%f", &options.snr) != 1) {
          fprintf(stderr, "%s: invalid sync SNR `%s'\n", argv[0], optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'c':
        if (tolower(*optarg) == 'i') {
          options.q_channel = false;
        } else if (tolower(*optarg) == 'q') {
          options.q_channel = true;
        } else {
          fprintf(stderr, "%s: invalid component `%c'\n", argv[0], *optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'g':
        if ((options.guess = strdup(optarg)) == NULL) {
          fprintf(stderr, "%s: memory exhausted\n", argv[0]);
          exit(EXIT_FAILURE);
        }
        break;

      case 'f':
        if ((options.file = strdup(optarg)) == NULL) {
          fprintf(stderr, "%s: memory exhausted\n", argv[0]);
          exit(EXIT_FAILURE);
        }
        break;
      
      case 'o':
        if ((options.ofile = strdup(optarg)) == NULL) {
          fprintf(stderr, "%s: memory exhausted\n", argv[0]);
          exit(EXIT_FAILURE);
        }
        break;
      
      case 'l':
        if (sscanf(optarg, "%u", &options.port) != 1 
            || options.port == 0
            || options.port > 65536) {
          fprintf(stderr, "%s: invalid port `%s'\n", argv[0], optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'i':
        if (sscanf(optarg, "%u", &options.iters) != 1) {
          fprintf(stderr, "%s: invalid iteration number `%s'\n", argv[0], optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'h':
        help(argv[0]);
        exit(EXIT_SUCCESS);

      case '?':
        exit(EXIT_FAILURE);
    }
  }

  if (!run(argv[0], &options))
    exit(EXIT_FAILURE);

  return 0;
}

# ccsds-tc
`ccsds-tc` is a small project that attempts to systematize the decoding of space packets as received by ground stations of the Amateur DSN. 
This is basically EB3FRN's fault for talking me into this :P

The project includes a custom implementation of small BCJR-based decoder for CCSDS turbo codes (according to CCSDS 131.0-B-3), a descrambler 
and the logic to compute the CRC16 of decoded frames.

## Build
This is the typical CMake project, just change to build, run `cmake` and `make`:

```
% cd build
% cmake ..
% make
```

## Running the project
The program you may be looking for is `ccsds-tool`. It supports 2 execution modes:
  * **Guess mode**, in which the tool will open an already demodulated baseband signal at the channel's symbol rate and attempt to blindly guess the channel parameters, and
  * **Decoding mode**, in which you pass the decoder parameters to tool in the command line, and it attempts to decode the input in real time.
  
By default, the tool operates in decoding mode, reading 32-bit complex float samples from the standard input and writing decoded (and CRC-checked) frames
to the standard output. Samples can be read from a file or a TCP socket as well:

```
% ./ccsds-tool -f samples-LLR.raw -o frames.bin
CCSDS tool v0.1 for the Amateur DSN by EA1IYR
(c) 2021 Gonzalo J. Carracedo - https://actinid.org
  Code rate:       1/6
  Block length:    8920 bits
  Turbocode iters: 1
  Channel:         Q
  Sync SNR:        19.03 dB
  Input file: /tmp/samples-LLR.raw
  Output file: frames.bin
Decode rate: 882.88 kbps
./ccsds-tool: 82510 bytes decoded
```

Guess mode operates from existing files only and it is activated by means of the command line option `-g`:

```
% ./ccsds-tool -g samples-LLR.raw 
Looking for syncwords for r = 1/2
Looking for syncwords for r = 1/3
Looking for syncwords for r = 1/4
Looking for syncwords for r = 1/6
Candidate CCSDS turbocode found
    Code rate:    1/6
    Channel Q:    94 frames
    Frame length: 8920 bit
```

Details on other tool options can be obtained by running `ccsds-tool --help`


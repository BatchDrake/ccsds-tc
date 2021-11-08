cmake_minimum_required(VERSION 3.10)

set(
  CCSDS_TC_SOURCES
  ccsds-tc.c)

set(
  CCSDS_TC_HEADERS
  ccsds-tc.h)

project(ccsds_tc VERSION 0.1)

add_executable(
  ccsds-tc-test
  ${CCSDS_TC_SOURCES} 
  ${CCSDS_TC_HEADERS} test.c)

add_executable(
  ccsds-tc-test-fast
  ${CCSDS_TC_SOURCES} 
  ${CCSDS_TC_HEADERS} test.c)

target_link_libraries(
  ccsds-tc-test
  m)

target_link_libraries(
  ccsds-tc-test-fast
  m)

target_compile_options(
  ccsds-tc-test-fast PRIVATE 
  -DCCSDS_TC_INT_ARITHMETICS -DCCSDS_TC_NO_MAXSTAR -O9)

add_executable(
  ccsds-tool
  ${CCSDS_TC_SOURCES} 
  ${CCSDS_TC_HEADERS}
  ccsds-guess.c
  ccsds-guess.h
  ccsds-modem.c
  ccsds-modem.h
  ccsds-tool.c
  correlator.c
  correlator.h
  defs.h
  lfsr.c
  lfsr.h)

target_link_libraries(
  ccsds-tool
  m)

target_compile_options(
  ccsds-tool PRIVATE 
  -DCCSDS_TC_NO_MAXSTAR -O9)

cmake_minimum_required (VERSION 2.6)
project(rtl_power_fftw)
include(FindPkgConfig)

option(FLOAT_FFT "Perform FFTs with single precision floats")

pkg_check_modules(TCLAP REQUIRED tclap)
pkg_check_modules(RTLSDR REQUIRED librtlsdr)

include_directories(
  ${TCLAP_INCLUDE_DIRS}
  ${RTLSDR_INCLUDE_DIRS}
)

if(FLOAT_FFT)
  pkg_check_modules(FFTW3F REQUIRED fftw3f)
  include_directories(${FFTW3F_INCLUDE_DIRS})
  add_definitions(-DFLOAT_FFT)
else()
  pkg_check_modules(FFTW3 REQUIRED fftw3)
  include_directories(${FFTW3_INCLUDE_DIRS})
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE)

add_executable(rtl_power_fftw rtl_power_fftw.cxx)

target_link_libraries(rtl_power_fftw
  ${TCLAP_LDFLAGS}
  ${RTLSDR_LDFLAGS}
  m
)

if(FLOAT_FFT)
  target_link_libraries(rtl_power_fftw ${FFTW3F_LDFLAGS})
else()
  target_link_libraries(rtl_power_fftw ${FFTW3_LDFLAGS})
endif()

INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/rtl_power_fftw
  DESTINATION bin)

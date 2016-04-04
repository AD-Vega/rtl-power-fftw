cmake_minimum_required (VERSION 2.6)
project(rtl_power_fftw)
include(FindPkgConfig)

pkg_check_modules(TCLAP REQUIRED tclap)
pkg_check_modules(RTLSDR REQUIRED librtlsdr)
pkg_check_modules(FFTW3F REQUIRED fftw3f)

include_directories(
  ${TCLAP_INCLUDE_DIRS}
  ${RTLSDR_INCLUDE_DIRS}
  ${FFTW3F_INCLUDE_DIRS}
)

# Check for 'long long int' to make TCLAP happy.
include(CheckTypeSize)
check_type_size("long long int" LONG_LONG BUILTIN_TYPES_ONLY)
if(HAVE_LONG_LONG)
  add_definitions(-DHAVE_LONG_LONG)
else()
  message(FATAL_ERROR
    "${PROJECT_NAME} requires the 'long long int' type to compile!")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE)

add_executable(rtl_power_fftw
  src/acquisition.cxx
  src/dispatcher.cxx
  src/device.cxx
  src/interrupts.cxx
  src/params.cxx
  src/fft_worker.cxx
  src/output.cxx
  src/rtl_power_fftw.cxx
)

target_link_libraries(rtl_power_fftw
  ${FFTW3F_LDFLAGS}
  ${TCLAP_LDFLAGS}
  ${RTLSDR_LDFLAGS}
  pthread
  m
)

INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/rtl_power_fftw
  DESTINATION bin)

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/doc/rtl_power_fftw.1
  DESTINATION share/man/man1)

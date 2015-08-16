cmake_minimum_required (VERSION 2.6)
project(rtl_power_fftw)
include(FindPkgConfig)
pkg_check_modules(PKGCONFS REQUIRED tclap)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fopenmp ${PKGCONFS_CFLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE)

add_executable(rtl_power_fftw rtl_power_fftw.cxx)
target_link_libraries(rtl_power_fftw ${PKGCONFS_LDFLAGS} rtlsdr fftw3 m)

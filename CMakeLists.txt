cmake_minimum_required (VERSION 2.6)
project(rtl_power_fftw)

add_executable(rtl_power_fftw rtl_power_fftw.cxx)
target_link_libraries(rtl_power_fftw rtlsdr fftw3 m)

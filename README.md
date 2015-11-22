# rtl\_power\_fftw

`rtl_power_fftw` is a program that obtains a power spectrum from RTL
devices using the FFTW library to do FFT.

It is inspired by the program `rtl_power` in `librtlsdr`.  However, the
said program has several deficiencies that limit its usage in
demanding environments, such as radio astronomy. An inspection of
`rtl_power` in hope of modifying it and obtaining better performance
resulted in the conclusion that it would be an unfeasible task. 
Measurements of FFT performance showed that the leading library
in the field of FFT - `fftw` - makes mincemeat of the routine used in
`rtl_power`, even on simple processors such as raspberryPi. Therefore,
the following requirements for a program to obtain power spectrum from
rtlsdr devices were set out:

  - the new program should use the `fftw` library
  - it will (of course) use `librtlsdr` to interface device
  - it should process spectra in a separate thread from data acquisition,
	to optimize observation time (continuous sampling)
  - it should be friendly to use and adhere to the UNIX philosophy of
	only doing one thing, but doing it well
  - the output of the program should be easy to further use with the 
	standard tools (like `gnuplot`).
  
The desire to have simple code to handle option parsing lead to the 
choice of TCLAP and therefore C++. This further meant that to implement 
things in a neat way, C++11 functionality snuck into the program and 
therefore a modern, C++11 enabled compiler is needed to 
compile `rtl_power_fftw`.

## Installation

To compile the program, cd into the directory where you have cloned the code
and do:

    mkdir build
    cd build
    cmake ..
    make

This should make the `rtl_power_fftw` binary in the build directory.
If you copy it into a directory in your `PATH`, you can call it from everywhere.
You can also do `make install` and by default it will be copied to `/usr/local/bin`.

## Documentation

The man page of `rtl_power_fftw` is available [here](doc/rtl_power_fftw.1.md).
It is maintained in Markdown format suitable for reading online, but it is
also converted into a regular man page that gets installed along with the
program.

## TODO:

  - ... (there's always something, isn't it?!)

Many thanks to Andrej Lajovic for cleaning up the C++ code, implementing a
better buffer handling system and relentlessly improving the program.

# NAME

rtl_power_fftw - obtain a power spectrum from RTL2832U devices

# SYNOPSIS

rtl_power_fftw [`OPTION` ...]


# DESCRIPTION

`rtl_power_fftw` reads the data stream from a RTL2832U device, performs the Fast Fourier transform on it, converts the frequency-domain data into a power spectrum and (usually, but not necessarily) averages a number of such spectra to gain a better signal-to-noise ratio. The program is capable of continuous real-time data acquisition and processing: it uses the rtl-sdr library to access the RTL2832U device and the FFTW library to do the FFT. Data acquisition and Fast Fourier transform are done in separate program threads for maximum efficiency. The resulting power spectrum is written to the standard output.

Scans of wide frequency ranges are supported: if the range exceeds the device bandwidth, `rtl_power_fftw` can perform the measurement in multiple steps, changing the central frequency each time so that the desired frequency range is completely covered. See **FREQUENCY SCANNING** below for details.

The program can be stopped gracefully by sending the SIGINT signal to it (pressing Ctrl+C). If only one such signal is received, the program will finish the current frequency scan, write out the data and terminate. This is useful for quitting continuous acquisition mode (see `--continue`) as it ensures that the current scan is finished before the program stops. If two SIGINTs are received, the current acquisition will stop as soon as possible and no attempt will be done to finish the current scan (the data will still be written out, but some spectra might be missing from the last frequency scan). The third SIGINT will signal the operating system to terminate the program immediately.

`rtl_power_fftw` accepts various options and switches that determine how the data from the RTL device is acquired and processed.


# OPTIONS

`-B <file|->`,  `--baseline <file|->`
:   Subtract baseline, read baseline data from a file or from standard input. See **BASELINE AND WINDOW FUNCTION DATA** below for format and further considerations.

`-b <bins in FFT spectrum>`,  `--bins <bins in FFT spectrum>`
:    Number of bins in the FFT spectrum (must be an even number).

`--buffers <buffers>`
:   Number of read buffers: don't touch unless running out of memory.

`-c`,  `--continue`
:   Repeat the same measurement endlessly. The spectra are written sequentially to the output and a header is appended before each measurement (as usual).

`-d <device index>`,  `--device <device index>`
:   RTL-SDR device index of the device used for the measurement.

`-e <duration>`,  `--session-duration <duration>`
:   Similar to `--continue`, but with limited duration. The program will repeat the same set of measurements until the time is up. The ongoing frequency scan is always completed so the actual running time might exceed the specified duration, especially if large frequency spans are used. For argument format, see the `--time` option.

`-f <Hz|Hz:Hz>`,  `--freq <Hz|Hz:Hz>`
:   Center frequency of the receiver or the frequency range to scan. A number can be followed by a `k`, `M` or `G` multiplier, meaning that the frequency is expressed in kilohertz, megahertz or gigahertz. Frequency range consists of lower and upper bound, separated with colon.

`-g <1/10th of dB>`,  `--gain <1/10th of dB>`
:   Receiver gain, expressed in tenths of a decibel (e.g., 100 means 10 dB).

`-H <number>`,  `--threads <number>`
:   The number of concurrent threads for FFT. The default is to use one thread per CPU core.

`-l`,  `--linear`
:   Output linear power values instead of logarithmic.

`-n <repeats>`,  `--repeats <repeats>`
:   Number of spectra to average (incompatible with `-t`).

`-o <percent>`,  `--overlap <percent>`
:   Define lower boundary for overlap when frequency hopping (otherwise meaningless).

`-O <file>`,  `--output <file>`
:   The output file. If this option is not given, the output is written to stdout.

`-p <ppm>`,  `--ppm <ppm>`
:   Correct for the oscillator error of RTL-SDR device. The correction should be given in ppm.

`-q`,  `--quiet`
:   Limit verbosity. Use once for reduced verbosity, twice for an almost complete silence.

`-r <Hz>`,  `--rate <Hz>`
:   Sample rate of the receiver in Hz.

`-s <bytes>`,  `--buffer-size <bytes>`
:   Size of the read buffers (leave it as it is unless you know what you are doing).

`-T`,  `--strict-time`
:   End measurement when the time set with `--time` option is up, regardless of the number of gathered samples.

`-t <duration>`,  `--time <duration>`
:   Integration time (incompatible with -n). This is an _effective integration time_; see **INTEGRATION TIME** below for more info (in short, the measurement might take *longer* than that). The argument can be specified in days, hours, minutes, and seconds (or a combination thereof) by appending `d`, `h`, `m` or `s` suffixes to the numbers. It is permissible to omit the `s` for seconds.

`-w <file|->`,  `--window <file|->`
:   Use a window function, read data from a file or from standard input. See **BASELINE AND WINDOW FUNCTION DATA** below for format and further considerations.

`--`,  `--ignore_rest`
:   Ignore the rest of the labeled arguments following this flag.

`--version`
:   Display version information and exit.

`-h`,  `--help`
:   Display usage information and exit.


# OUTPUT FORMAT

Every spectrum that `rtl_power_fftw` writes to stdout is preceded by a few lines of metadata: these lines begin with a `#` character. Timestamps of data acquisition start and end are written, followed by a header line describing the data. This is an example output:

    # Acquisition start: 2015-11-22 17:59:34 UTC
    # Acquisition end: 2015-11-22 17:59:34 UTC
    #
    # frequency [Hz] power spectral density [dB/Hz]
    1.41940575e+09 -68.7714
    1.41940966e+09 -68.668
    ...

Data lines contain two columns: the frequency and the power spectral density. Note: the power values are expressed in logarithmic (dB) scale, but the reference power is not explicitly given. In particular, the numbers are **NOT** in dBm; they might be different for a different device and might also change (shift by a constant) in future versions of the program. In other words: don't assume anything. If you need absolute units, you have to calibrate your device against a known reference signal.

If several measurements are to be done, the consecutive spectra will be divided by blank lines (see **FREQUENCY SCANNING** for details).


# INTEGRATION TIME

The integration time specified with the `--time` option is usually considered to be the *effective* integration time, i.e., the total number of samples to be acquired divided by the sample rate. If all goes well, the program will run for about this long (plus some overhead for setting up the device etc.). However, if samples are dropped for any reason (for example, if the CPU can't perform the FFTs quickly enough to cope with the actual data rate), the program will run for whatever time required to collect the needed number of samples, which can be considerably longer than the effective integration time.

If you need the program to stop after a fixed time -- regardless of the actual number of samples collected -- use the `--strict-time` switch. Be warned, though, that only *acquisition* will be stopped after this time and it can take several more seconds for the FFT of the remaining data to be performed (this time overhead depends on the number of buffers used, see **BUFFERING** below).


# FREQUENCY SCANNING

If the frequency span is too large to be contained within a single measurement (i.e., it exceeds the device bandwidth), `rtl_power_fftw` will divide it into several consecutive measurements.

Of course, this raises a question: how to go about fitting several fixed-width (one device bandwidth) measurements into an arbitrary range? One could go for non-overlapping measurements, which yields data that is monotonously increasing in frequency, but then the whole scan might need to start *below* the lowest requested frequency, or end *above* the highest requested frequency, or even both. Even worse, these extended ranges could happen to contain frequencies not accepted by the device. Another approach is therefore used, namely to cover the requested frequency range exactly, but with overlapping measurements. Note that `rtl_power_fftw` will not make any presumptions on what to do with the overlaps: the overlapping spectra are simply written to the output and all further data treatment is up to the user. In case that your particular data treatment requires a certain minimum amount of overlap, you can use the option `--overlap` to set the desired lower bound for overlap in percentage of bandwidth.

All spectra within one scan of the desired frequency range are separated in the output by a single blank line. After the whole frequency range has been scanned, an additional blank line is printed, so the measurement *sets* are separated by two blank lines in total. This output format is directly suitable as an input for `gnuplot`.


# BASELINE AND WINDOW FUNCTION DATA

The expected input format for baseline and window function data is one value per line. If a line contains multiple values, the last (rightmost) value is used: this ensures that `rtl_power_fftw` can use its own output data as an input for baseline correction -- the frequency column is simply discarded. Lines starting with `#` are treated as comments and are ignored completely.

If both the baseline and window function data are to be read from standard input, the baseline data is read first, followed by the window function data.

The program does not check the window function data in any way, apart from the requirement to have precisely enough data points. Window function is only read in single precision, due to FFT being done with floats, and there is no need to overcomplicate things. Single precision FFT is faster than double on at least some hardware and more than precise enough, as input data is actually only 8-bit. Baseline data is in double precision, otherwise it would limit the precision of averaging arbitrarily huge number of spectra.


# BUFFERING

Upon starting, the program allocates several data buffers (five by default). At any given time, one of the buffers is used to store the incoming data from the device. When the buffer fills up, it is queued for processing by the FFT routine and an empty buffer is immediately taken to continue the data acquisition; at this point, the number of empty buffers is also recorded for statistical purposes (see below). If no buffers are empty, the data acquisition blocks until one of the buffers becomes available again. This is, of course, an unwanted scenario because it leads to dropped data.

At the end of the measurement, the program outputs a line with the statistics on the number of available (empty) buffers. This is an example of such a line:

    Buffer queue histogram: 0 0 0 6 34 1

The numbers report how many times a particular number of available buffers was encountered. The first number corresponds to zero available buffers, the next one to one available buffer and so on. In this particular case, at least three buffers were available at all times: three buffers were available on six occurrences, four buffers were available on 34 occurrences and all five buffers were only available once (when the program started and there was no data yet).

As long as the first number remains zero, you are fine - there was no data loss. If the first number happens to be nonzero and also exceeds the other numbers, this means that your CPU is too slow and cannot perform the FFTs quickly enough to match the incoming data rate. You might be better off with a smaller FFT size or a slower sampling rate.

On the contrary, if the first number is nonzero but is relatively small compared to the other numbers, it might simply mean that the available CPU power fluctuates heavily (e.g., if you have a fast processor but other CPU-intensive tasks are running at the same time). In such a case, you can try increasing the number of buffers with the `--buffers` option and see if that helps.

Another scenario occurs if you have enough computing power but the memory is limited: in such a (rare) case, you might actually want to *reduce* the number of buffers.

The size (length) of the buffers is computed automatically to best match the requirements of the measurement. This is the recommended practice in most circumstances. However, if you feel that you have a very good reason to fiddle with the buffer size, you can do so with the `--buffer-size` option. But do keep in mind that the buffer size should be a multiple of 16384 (this is a requirement of the rtl-sdr library).


# EXAMPLES

A basic call to `rtl_power_fftw` might look like this:

    rtl_power_fftw -f 1420405752 -t 10 -b 512 > spectrum.dat

This will set the central frequency of the receiver to 1420405752 Hz (the frequency of the hydrogen line), use a 512-point FFT to transform the acquired signal, average the data for ten seconds and dump the averaged spectrum to a file named _spectrum.dat_.

By the virtue of the output data being suitable for direct use in `gnuplot`, the following pipeline can be used to acquire a spectrum and draw it into a PNG image (for variety, the `-n` option is used this time to request the average of 100 spectra):

    rtl_power_fftw -f 1420405752 -n 100 -b 512 |\
       gnuplot -e "set term png; unset key; plot '-' w l" >plot.png

For quick-and-dirty live monitoring, you can do:

    rtl_power_fftw -f 1420405752 -n 100 -b 512 -c |\
       sed -u '/rtl-power-fftw/s/.*/plot "-"/;/^$/{N;s/^\n$/e/}' |\
       gnuplot

In this pipeline, `sed` intervenes by replacing the header and separators written by `rtl_power_fftw` with inline commands for `gnuplot`.

To scan frequencies between 100 MHz and 110 MHz and subtract baseline data from each scan, you could do:

    rtl_power_fftw -f 100M:110M -B baseline_data.dat > spectrum.dat
    
This example also illustrates the fact that for all the options where it is possible, the program selects some safe default values and the options can be omitted. Although be noted that omiting the option to specify number of bins (`-b`) and relying on its default value while subtracting baseline is a discouraged practise. You should always specify `--bins` along with `--baseline`.


## AUTHORS

Klemen Blokar <klemen.blokar@ad-vega.si>  
Andrej Lajovic <andrej.lajovic@ad-vega.si>

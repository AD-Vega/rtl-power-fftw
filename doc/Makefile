# This recipe generates the man page of rtl_power_fftw from the Markdown
# source. It is intended to be run manually when the Markdown source changes
# and the resulting roff file is then distributed along with the program.
# This way, people who just want to build and install rtl_power_fftw are not
# required to have pandoc installed as well (as it is quite a hefty
# package).

MAN_FILE = rtl_power_fftw.1
MD_SOURCE = $(MAN_FILE).md

$(MAN_FILE): $(MD_SOURCE)
	pandoc -s \
	       -t man \
	       -V title=rtl_power_fftw \
	       -V section=1 \
	       $(MD_SOURCE) \
	       -o $(MAN_FILE)

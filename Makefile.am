AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  Chi2checker.h

ROOTDICTS = 

pcmdir = $(libdir)
nobase_dist_pcm_DATA = 

lib_LTLIBRARIES = \
  libchi2checker_io.la \
  libchi2checker.la

libchi2checker_io_la_SOURCES = \
  $(ROOTDICTS) 

libchi2checker_la_SOURCES = \
  Chi2checker.cc

libchi2checker_io_la_LIBADD = \
  -lphool -lSubsysReco -lmbd_io -ljetbase

libchi2checker_la_LIBADD = \
  libchi2checker_io.la \
  -lcalo_io \
  -lSubsysReco -lphg4hit -lffaobjects -lmbd_io \
  -ljetbase -lglobalvertex -lcentrality -lcalotrigger -lphparameter
#-lg4jets_io
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_io \
  testexternals

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD   = libchi2checker_io.la

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libchi2checker.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)

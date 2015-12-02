PYTHIA8 = $(PYTHIA)

-include $(PYTHIA8)/share/Pythia8/examples/Makefile.inc
ifeq ($(GZIP_USE),true)
  CXX_COMMON+= -DGZIPSUPPORT
  CXX_COMMON+= -L$(BOOST_LIB) -Wl,-rpath $(BOOST_LIB) -lboost_iostreams
  CXX_COMMON+= -L$(GZIP_LIB) -Wl,-rpath $(GZIP_LIB) -lz
endif
ifneq ("$(wildcard $(PYTHIA8)/lib/libpythia8.a)","")
  PREFIX_LIB=$(PYTHIA8)/lib
  PREFIX_INCLUDE=$(PYTHIA8)/include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON) -Wl,-rpath $(PREFIX_LIB) -ldl
#CXX_COMMON:=-I$(PREFIX_INCLUDE) $(PREFIX_LIB)

.SECONDEXPANSION:
.PHONY: all clean

Makefile.inc:
	$(error Error: PYTHIA must be configured, please run "./configure"\
                in the top PYTHIA directory)

$(PREFIX_LIB)/libpythia8.a :
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)

Evis : Evis.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON)

hadron_dist : hadron_dist.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON)

hadron_dist_moroi : hadron_dist_moroi.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON)

hadron_dist_lhe : hadron_dist_lhe.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON)

hadron_Evis_lhe : hadron_Evis_lhe.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON)

clean:
	@rm -f main[0-9][0-9]; rm -f out[0-9][0-9]; rm -f weakbosons.lhe;\
	rm -f mymain[0-9][0-9]; rm -f myout[0-9][0-9]; rm -f hist.root;\
	rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.*

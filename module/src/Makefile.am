SUBDIRS = distributions

jagsmod_LTLIBRARIES = ZIBGeometric.la

ZIBGeometric_la_SOURCES = ZIBGeometric.cc

ZIBGeometric_la_CPPFLAGS = -I$(includedir)/JAGS 

ZIBGeometric_la_LIBADD = distributions/ZIBGeometricdist.la
if WINDOWS
ZIBGeometric_la_LIBADD += -ljags-$(JAGS_MAJOR)
else
ZIBGeometric_la_LIBADD += -ljags
endif

ZIBGeometric_la_LDFLAGS = -module -avoid-version
if WINDOWS
ZIBGeometric_la_LDFLAGS += no-undefined
endif
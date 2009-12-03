#-----------------------------------------------------------------------------
# GNU makefile for building WCSLIB 4.4.
#
# Summary of the main targets
# ---------------------------
#   all:       Do 'make all' in each subdirectory.
#   check:     Do 'make check' in each subdirectory.
#   tests:     Do 'make tests' in each subdirectory.
#   install:   Do 'make install' in each subdirectory.
#   realclean (or cleanest): Do 'make realclean' in each subdirectory.
#   show_all:  Print the values of all variables used in this GNUmakefile.
#   writable:  Run chmod recursively to make all sources writable.
#
# Notes:
#   1) If you need to make changes then preferably modify makedefs.in instead.
#
#   2) Refer also to the GNUmakefiles in subdirectories, particularly
#      C/GNUmakefile.
#
# Author: Mark Calabretta, Australia Telescope National Facility
# http://www.atnf.csiro.au/~mcalabre/index.html
# $Id: GNUmakefile,v 4.4.1.6 2009/09/14 06:49:26 cal103 Exp cal103 $
#-----------------------------------------------------------------------------
# Get configure settings.
include makedefs

.PHONY : build check chmod cleanest install realclean show_all tests writable

build :
	-@ for DIR in $(SUBDIRS) ; do \
	     echo '' ; \
	     $(TIMER) ; \
	     $(MAKE) -k -C $$DIR build ; \
	   done

check tests : show_all
	-@ echo ''
	-@ $(TIMER)
	 @ for DIR in $(SUBDIRS) ; do \
	     echo '' ; \
	     $(MAKE) -i -C $$DIR cleaner ; \
	   done
	-@ echo ''
	 @ for DIR in $(TSTDIRS) ; do \
	     echo '' ; \
	     $(TIMER) ; \
	     $(MAKE) -k -C $$DIR $@ ; \
	   done

install :
	 @ for DIR in $(INSTDIR) ; do \
	     $(MAKE) -k -C $$DIR $@ ; \
	   done
	   $(INSTALL) -m 444 wcsconfig.h wcsconfig_f77.h $(INCDIR)
	   $(INSTALL) -m 444 CHANGES $(INCDIR)/CHANGES
	   $(INSTALL) -m 444 COPYING $(INCDIR)/COPYING
	   $(INSTALL) -m 444 README  $(INCDIR)/README
	   if [ ! -d "$(LIBDIR)/pkgconfig" ] ; then \
	     $(INSTALL) -d -m 2775 $(LIBDIR)/pkgconfig ; \
	   fi
	   $(INSTALL) -m 444 wcslib.pc $(LIBDIR)/pkgconfig/wcslib.pc

realclean cleanest :
	for DIR in $(SUBDIRS) doxygen ; do \
	   $(MAKE) -C $$DIR realclean ; \
	done
	$(RM) -r autom4te.cache autoscan.log
	$(RM) confdefs.h conftest.*
	$(RM) config.log config.status configure.lineno
	$(RM) makedefs wcslib.pc
	$(RM) wcsconfig.h wcsconfig_*.h
	$(RM) wcslib-*.tar.gz

show_all ::
	-@ echo 'Subdirectories to be built...'
	-@ echo '  SUBDIRS     := $(SUBDIRS)'
	-@ echo '  TSTDIRS     := $(TSTDIRS)'
	-@ echo ''

writable :
	  chmod -R u+w .

GNUmakefile : makedefs ;

makedefs : makedefs.in config.status
	-@ echo ''
	-@ $(TIMER)
	   ./config.status

config.status : configure
	-@ echo ''
	-@ $(TIMER)
	-@ echo ''
	-@ echo "Environment variables that affect 'configure':"
	-@ echo "  CC       = $${CC-(undefined)}"
	-@ echo "  CFLAGS   = $${CFLAGS-(undefined)}"
	-@ echo "  CPP      = $${CPP-(undefined)}"
	-@ echo "  CPPFLAGS = $${CPPFLAGS-(undefined)}"
	-@ echo "  F77      = $${F77-(undefined)}"
	-@ echo "  FFLAGS   = $${FFLAGS-(undefined)}"
	-@ echo "  LDFLAGS  = $${LDFLAGS-(undefined)}"
	-@ echo ''
	   ./configure --no-create


#-----------------------------------------------------------------------------
# These are for code management.

.PHONY : dist

dist :
	   $(MAKE) -C doxygen cleanest build
	   $(MAKE) realclean
	-@ echo $(WCSLIBV)/C/RCS        >  wcslib.X
	-@ echo $(WCSLIBV)/C/flexed/RCS >> wcslib.X
	-@ echo $(WCSLIBV)/C/test/RCS   >> wcslib.X
	-@ echo $(WCSLIBV)/doxygen/RCS  >> wcslib.X
	-@ echo $(WCSLIBV)/Fortran/RCS  >> wcslib.X
	-@ echo $(WCSLIBV)/Fortran/test/RCS >> wcslib.X
	-@ echo $(WCSLIBV)/makedefs     >> wcslib.X
	-@ echo $(WCSLIBV)/other        >> wcslib.X
	-@ echo $(WCSLIBV)/pgsbox/RCS   >> wcslib.X
	-@ echo $(WCSLIBV)/releases     >> wcslib.X
	-@ echo $(WCSLIBV)/RCS          >> wcslib.X
	-@ echo $(WCSLIBV)/TODO         >> wcslib.X
	-@ echo $(WCSLIBV)/utils/RCS    >> wcslib.X
	-@ echo $(WCSLIBV)/wcslib.T     >> wcslib.X
	-@ echo $(WCSLIBV)/wcslib.X     >> wcslib.X
	   rm -f $(WCSLIBV).tar.bz2
	   tar cf - -C .. -X wcslib.X $(WCSLIBV) | \
	     tar t | \
	     grep -v '/$$' | \
	     sort > wcslib.T
	   rm -f wcslib.X
	   tar cvf $(WCSLIBV).tar -C .. -T wcslib.T
	   rm -f wcslib.T
	   bzip2 $(WCSLIBV).tar
	   chmod 444 $(WCSLIBV).tar.bz2

install_dist :
	   mv -f $(WCSLIBV).tar.bz2 ./releases/
	   cp -fp ./releases/$(WCSLIBV).tar.bz2 /nfs/ftp/software/wcslib/
	   (cd /nfs/ftp/software/wcslib/ && \
	     rm -f wcslib.tar.bz2 && \
	     ln -s $(WCSLIBV).tar.bz2 wcslib.tar.bz2)
	   cp -fp CHANGES wcslib.pdf ~/public_html/WCS/
	   rsync --archive --delete html/ ~/public_html/WCS/wcslib/

configure : configure.ac
	-@ echo ''
	-@ $(TIMER)
	   autoconf

# Code development overrides must be included specifically before 'configure'
# generates makedefs.
-include flavours

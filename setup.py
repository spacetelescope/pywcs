#!/usr/bin/env python

CONTACT = "mdroe@stsci.edu"

import sys
from distutils.core import setup, Extension
from os.path import join

######################################################################
# NUMPY
try:
    import numpy
except ImportError:
    print "numpy must be installed to build pywcs."
    print "ABORTING."
    sys.exit(1)

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

######################################################################
# PyFITS
try:
    import pyfits
except ImportError:
    print "WARNING: PyFITS must be installed to use pywcs."
    print "         Since this is not a build-time dependency, the build will proceed."

######################################################################
# WCSLIB
WCSVERSION = "4.3"
WCSLIB = "wcslib-%s" % WCSVERSION # Path to wcslib
WCSLIBC = join(WCSLIB, "C") # Path to wcslib source files
WCSFILES = [ # List of wcslib files to compile
    'flexed/wcsulex.c',
    'flexed/wcsutrn.c',
    'cel.c',
    'lin.c',
    'log.c',
    'prj.c',
    'spc.c',
    'sph.c',
    'spx.c',
    'tab.c',
    'wcs.c',
    'wcsfix.c',
    'wcshdr.c',
    'wcspih.c',
    'wcstrig.c',
    'wcsunits.c',
    'wcsutil.c']
WCSFILES = [join(WCSLIBC, x) for x in WCSFILES]

######################################################################
# WCSLIB CONFIGURATION

# The only configuration parameter needed at compile-time is how to
# specify a 64-bit signed integer.  Python's ctypes module can get us
# that information, but it is only available in Python 2.5 or later.
# If we can't be absolutely certain, we default to "long long int",
# which is correct on most platforms (x86, x86_64).  If we find
# platforms where this heuristic doesn't work, we may need to hardcode
# for them.
def determine_64_bit_int():
    try:
        import ctypes
    except ImportError:
        print "WARNING: Unable to determine a suitable 64-bit integer type."
        print "         Defaulting to 'long long int', which should work with most"
        print "         platforms, but your build may be broken."
        print "         Please contact <%s> with details about your platform." % CONTACT
        return "long long int"

    if ctypes.sizeof(ctypes.c_longlong) == 8:
        return "long long int"
    elif ctypes.sizeof(ctypes.c_long) == 8:
        return "long int"
    elif ctypes.sizeof(ctypes.c_int) == 8:
        return "int"
    else:
        print "WARNING: Could not find a suitable 64-bit integer type."
        print "         Defaulting to 'long long int', but your build may be broken."
        print "         Please contact <%s> with details about your platform." % CONTACT
        return "long long int"

fd = open("src/wcsconfig.h", "w")
fd.write("""
/* WCSLIB library version number. */
#define WCSLIB_VERSION %s

/* 64-bit integer data type. */
#define WCSLIB_INT64 %s
""" % (WCSVERSION, determine_64_bit_int()))
fd.close()

######################################################################
# GENERATE DOCSTRINGS IN C
docstrings = {}
execfile("doc/docstrings.py", docstrings)
keys = docstrings.keys()
keys.sort()
fd = open("src/docstrings.h", "w")
fd.write('/* This file is autogenerated by setup.py.  To edit its contents\n')
fd.write('   edit doc/docstrings.py\n')
fd.write('*/\n')
for key in keys:
    if key.startswith('__'):
        continue
    val = docstrings[key].lstrip().encode("string_escape").replace('"', '\\"')
    fd.write("static const char doc_%s[] = \"%s\";\n\n" % (key, val))
fd.close()

######################################################################
# DISTUTILS SETUP
setup(name="pywcs",
      version="1.0a1-%s" % WCSVERSION,
      description="Python wrappers to WCSLIB",
      author="Michael Droettboom",
      author_email=CONTACT,
      url="http://projects.scipy.org/astropy/astrolib/wiki/WikiStart",
      packages=['pywcs'],
      ext_modules=[
        Extension('pywcs._pywcs',
                  ['src/pywcs.c'] +
                  WCSFILES,
                  include_dirs=[
                    numpy_include,
                    WCSLIBC,
                    "src"
                    ]
                  )
      ]
      )

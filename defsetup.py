#!/usr/bin/env python

CONTACT = "Michael Droettboom"
EMAIL = "mdroe@stsci.edu"

import sys
from distutils.core import setup, Extension
from os.path import join
import os.path

######################################################################
# CONFIGURATION
DEBUG = False
OPENMP = False

######################################################################
# NUMPY
try:
    import numpy
except ImportError:
    print "numpy must be installed to build pywcs."
    print "ABORTING."
    raise

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
    print "         Since this is not a build-time dependency, the "
    print "         build will proceed."

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
        try:
            import ctypes
        except ImportError:
            raise ValueError()

        if ctypes.sizeof(ctypes.c_longlong) == 8:
            return "long long int"
        elif ctypes.sizeof(ctypes.c_long) == 8:
            return "long int"
        elif ctypes.sizeof(ctypes.c_int) == 8:
            return "int"
        else:
            raise ValueError()

    except ValueError:
        print "WARNING: Could not find a suitable 64-bit integer type."
        print "         Defaulting to 'long long int', but your build may be broken."
        print "         Please contact <%s> with details about your platform." % EMAIL
        return "long long int"

if os.path.exists("pywcs"):
    srcroot = 'pywcs'
else:
    srcroot = '.'
fd = open(join(srcroot, 'src', 'wcsconfig.h'), "w")
fd.write("""
/* WCSLIB library version number. */
#define WCSLIB_VERSION %s

/* 64-bit integer data type. */
#define WCSLIB_INT64 %s
""" % (WCSVERSION, determine_64_bit_int()))
fd.close()

######################################################################
# GENERATE DOCSTRINGS IN C
sys.path.append(join('.', srcroot, "lib"))
docstrings = {}
execfile(join(srcroot, 'doc', 'docstrings.py'), docstrings)
keys = docstrings.keys()
keys.sort()
fd = open(join(srcroot, 'src', 'docstrings.h'), "w")
fd.write("""/*
DO NOT EDIT!

This file is autogenerated by setup.py.  To edit its contents,
edit doc/docstrings.py
*/

""")
for key in keys:
    if key.startswith('__'):
        continue
    val = docstrings[key].lstrip().encode("string_escape").replace('"', '\\"')
    fd.write('/*@unused@*/ static const char doc_%s[] = "%s";\n\n' % (key, val))
fd.close()

######################################################################
# PYWCS-SPECIFIC AND WRAPPER SOURCE FILES
PYWCS_SOURCES = [ # List of pywcs files to compile
    'distortion.c',
    'distortion_wrap.c',
    'pipeline.c',
    'pyutil.c',
    'pywcs.c',
    'sip.c',
    'sip_wrap.c',
    'str_list_proxy.c',
    'wcslib_wrap.c']
PYWCS_SOURCES = [join('src', x) for x in PYWCS_SOURCES]

######################################################################
# DISTUTILS SETUP
libraries = []
define_macros = [('ECHO', None)]
undef_macros = []
extra_compile_args = []
if DEBUG:
    define_macros.append(('DEBUG', None))
    undef_macros.append('NDEBUG')
    if not sys.platform.startswith('sun'):
        extra_compile_args.extend(["-fno-inline", "-O0", "-g"])
else:
    # Define ECHO as nothing to prevent spurious newlines from
    # printing within the libwcs parser
    define_macros.append(('NDEBUG', None))
    undef_macros.append('DEBUG')

if not sys.platform.startswith('sun'):
    if OPENMP:
        extra_compile_args.append('-fopenmp')
        libraries.append('gomp')
    else:
        extra_compile_args.append('-Wno-unknown-pragmas')

PYWCS_EXTENSIONS = [Extension('pywcs._pywcs',
                  WCSFILES + PYWCS_SOURCES,
                  include_dirs=[
                    numpy_include,
                    join(srcroot, WCSLIBC),
                    WCSLIBC,
                    join(srcroot, "src")
                    ],
                  define_macros=define_macros,
                  undef_macros=undef_macros,
                  extra_compile_args=extra_compile_args,
                  libraries=libraries
                  )
        ]

pkg = ["pywcs", "pywcs.include", "pywcs.include.wcslib"]

setupargs = {
    'version' :	    "1.4.1-%s" % WCSVERSION,
    'description':  "Python wrappers to WCSLIB",
    'author' :      CONTACT,
    'author_email': EMAIL,
    'url' :         "http://projects.scipy.org/astropy/astrolib/wiki/WikiStart",
    'platforms' :			["unix","windows"],
    'ext_modules' :			PYWCS_EXTENSIONS,
    'package_dir' : {
        'pywcs': 'lib',
        'pywcs.include': 'src',
        'pywcs.include.wcslib': WCSLIBC},
    'package_data' : {
        'pywcs.include': ['*.h'],
        'pywcs.include.wcslib': ['*.h']}
}


#!/usr/bin/env python

from __future__ import division # confidence high

CONTACT = "Michael Droettboom"
EMAIL = "mdroe@stsci.edu"

import sys
from distutils.core import setup, Extension
from os.path import join
import os.path

######################################################################
# CONFIGURATION
# BUILD may be 'debug', 'profile', or 'release'
BUILD = 'release'
OPENMP = False

######################################################################
# NUMPY
try:
    import numpy
except ImportError:
    print "numpy must be installed to build pywcs."
    print "ABORTING."
    raise

major, minor, rest = numpy.__version__.split(".", 2)
if (major, minor) < (1, 3):
    print "numpy version 1.3 or later must be installed to build pywcs."
    print "ABORTING."
    raise ImportError

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

######################################################################
# WCSLIB
WCSVERSION = "4.3.3"
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
def encode_docstring(s):
    return s.lstrip().encode("string_escape").replace('"', '\\"')

sys.path.append(join('.', srcroot, "lib"))
docstrings = {}
execfile(join(srcroot, 'doc', 'docstrings.py'), docstrings)
keys = [key for key in docstrings.keys() if not key.startswith('__')]
keys.sort()
for key in keys:
    docstrings[key] = docstrings[key].lstrip() + '\0'
fd = open(join(srcroot, 'src', 'docstrings.h'), "w")
fd.write("""/*
DO NOT EDIT!

This file is autogenerated by setup.py.  To edit its contents,
edit doc/docstrings.py
*/

#ifndef __DOCSTRINGS_H__
#define __DOCSTRINGS_H__

void fill_docstrings(void);

""")
for key in keys:
    val = docstrings[key]
    fd.write('extern char doc_%s[%d];\n' % (key, len(val)))
fd.write("\n#endif\n\n")
fd.close()

fd = open(join(srcroot, 'src', 'docstrings.c'), "w")
fd.write("""/*
DO NOT EDIT!

This file is autogenerated by setup.py.  To edit its contents,
edit doc/docstrings.py
*/

#include <string.h>
#include "docstrings.h"

""")
for key in keys:
    val = docstrings[key]
    fd.write('char doc_%s[%d];\n' % (key, len(val)))

fd.write("\nvoid fill_docstrings(void)\n{\n")
for key in keys:
    val = docstrings[key]
    # For portability across various compilers, we need to fill the
    # docstrings in 256-character chunks
    for i in range(0, len(val), 256):
        chunk = val[i:i+256].encode("string_escape").replace('"', '\\"')
        fd.write('strncpy(doc_%s + %d, "%s", %d);\n' % (
                key, i, chunk, min(len(val) - i, 256)))
    fd.write("\n")
fd.write("\n}\n\n")
fd.close()


######################################################################
# PYWCS-SPECIFIC AND WRAPPER SOURCE FILES
PYWCS_VERSION = '1.6'
VERSION = '%s-%s' % (PYWCS_VERSION, WCSVERSION)
PYWCS_SOURCES = [ # List of pywcs files to compile
    'distortion.c',
    'distortion_wrap.c',
    'docstrings.c',
    'pipeline.c',
    'pyutil.c',
    'pywcs.c',
    'pywcs_api.c',
    'sip.c',
    'sip_wrap.c',
    'str_list_proxy.c',
    'wcslib_wrap.c']
PYWCS_SOURCES = [join('src', x) for x in PYWCS_SOURCES]

######################################################################
# DISTUTILS SETUP
libraries = []
define_macros = [('ECHO', None), ('WCSTRIG_MACRO', None),
                 ('PYWCS_BUILD', None), ('_GNU_SOURCE', None)]
undef_macros = []
extra_compile_args = []
if BUILD.lower() == 'debug':
    define_macros.append(('DEBUG', None))
    undef_macros.append('NDEBUG')
    if not sys.platform.startswith('sun') and \
       not sys.platform == 'win32':
        extra_compile_args.extend(["-fno-inline", "-O0", "-g"])
elif BUILD.lower() == 'profile':
    define_macros.append(('NDEBUG', None))
    undef_macros.append('DEBUG')
    if not sys.platform.startswith('sun'):
        extra_compile_args.extend(["-O3", "-g"])
elif BUILD.lower() == 'release':
    # Define ECHO as nothing to prevent spurious newlines from
    # printing within the libwcs parser
    define_macros.append(('NDEBUG', None))
    undef_macros.append('DEBUG')
else:
    raise ValueError("BUILD should be one of 'debug', 'profile', or 'release'")

if sys.platform == 'win32':
    define_macros.append(('YY_NO_UNISTD_H', None))
    define_macros.append(('_CRT_SECURE_NO_WARNINGS', None))

if not sys.platform.startswith('sun') and \
   not sys.platform == 'win32':
    if OPENMP:
        extra_compile_args.append('-fopenmp')
        libraries.append('gomp')
    else:
        extra_compile_args.extend(['-Wno-unknown-pragmas'])

PYWCS_EXTENSIONS = [
    Extension('pywcs._pywcs',
              WCSFILES + PYWCS_SOURCES,
              include_dirs =
              [numpy_include,
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

pkg = ["pywcs" ]

setupargs = {
    'version' :	    VERSION,
    'description':  "Python wrappers to WCSLIB",
    'author' :      CONTACT,
    'author_email': EMAIL,
    'url' :         "http://projects.scipy.org/astropy/astrolib/wiki/WikiStart",
    'platforms' :			["unix","windows"],
    'ext_modules' :			PYWCS_EXTENSIONS,
    'data_files' : [
                    ( 'pywcs/include', ['src/*.h']),
                    ( 'pywcs/include/wcslib', [ WCSLIBC + '/*.h'] ),
                    ],
}


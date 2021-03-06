#!/usr/bin/env python

from __future__ import with_statement, division # confidence high

CONTACT = "Michael Droettboom"
EMAIL = "mdroe@stsci.edu"

from distutils.core import setup, Extension
from distutils.dist import Distribution
from distutils.errors import DistutilsError
import glob
from os.path import join
import os.path
import shutil
import sys

if os.path.exists("pywcs"):
    srcroot = 'pywcs'
else:
    srcroot = '.'
sys.path.append(join('.', srcroot, "lib/pywcs"))
sys.path.append('.')

def b(s):
    return s.encode('ascii')

if sys.version_info[0] >= 3:
    def string_escape(s):
        s = s.decode('ascii').encode('ascii', 'backslashreplace')
        s = s.replace(b('\n'), b('\\n'))
        return s.decode('ascii')
    from io import StringIO
    string_types = (str, bytes)
else:
    def string_escape(s):
        return s.encode('string_escape')
    from cStringIO import StringIO
    string_types = (str, unicode)

######################################################################
# CONFIGURATION
# BUILD may be 'debug', 'profile', or 'release'
BUILD = 'release'
OPENMP = False

######################################################################
# Helper class
def write_if_different(filename, data):
    data = data.encode('ascii')

    if os.path.exists(filename):
        with open(filename, 'rb') as fd:
            original_data = fd.read()
    else:
        original_data = None

    if original_data != data:
        with open(filename, 'wb') as fd:
            fd.write(data)

######################################################################
# NUMPY
try:
    import numpy
except ImportError:
    print("numpy must be installed to build pywcs.")
    print("ABORTING.")
    raise

major, minor, rest = numpy.__version__.split(".", 2)
if (int(major), int(minor)) < (1, 3):
    print("numpy version 1.3 or later must be installed to build pywcs.")
    print("ABORTING.")
    raise ImportError

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

######################################################################
# WCSLIB
WCSVERSION = "4.10"
WCSLIB = "wcslib" # Path to wcslib
WCSLIB_PATCHED = "wcslib"
WCSLIBC = join(WCSLIB_PATCHED, "C") # Path to wcslib source files
WCSFILES = [ # List of wcslib files to compile
    'flexed/wcsbth.c',
    'flexed/wcspih.c',
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
    'wcserr.c',
    'wcsfix.c',
    'wcshdr.c',
    'wcsprintf.c',
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

h_file = StringIO()
h_file.write("""
/* WCSLIB library version number. */
#define WCSLIB_VERSION %s

/* 64-bit integer data type. */
#define WCSLIB_INT64 %s

/* Windows needs some other defines */
/* I think the only one we need is _WIN32, but the others don't hurt - Mark S. 2012-02-14 */
#if defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__) || defined (__MINGW64__)

/* */
#ifndef YY_NO_UNISTD_H
#define YY_NO_UNISTD_H
#endif

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifndef _NO_OLDNAMES    /* for mingw32 */
#define _NO_OLDNAMES
#endif

#ifndef NO_OLDNAMES     /* for mingw64 */
#define NO_OLDNAMES
#endif

#ifndef __STDC__        /* for MS Visual C */
#define __STDC__ 1
#endif

#endif
""" % (WCSVERSION, determine_64_bit_int()))
write_if_different(join(srcroot, 'src', 'wcsconfig.h'), h_file.getvalue())

######################################################################
# GENERATE DOCSTRINGS IN C
docstrings = {}
with open(join(srcroot, 'doc', 'docstrings.py'), 'rb') as fd:
    docstrings_content = fd.read()
exec(docstrings_content, docstrings)
keys = [key for key in docstrings.keys()
        if not key.startswith('__') and type(key) in string_types]
keys.sort()
for key in keys:
    docstrings[key] = docstrings[key].encode('utf8').lstrip()

h_file = StringIO()
h_file.write("""/*
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
    h_file.write('extern char doc_%s[%d];\n' % (key, len(val)))
h_file.write("\n#endif\n\n")

write_if_different(join(srcroot, 'src', 'docstrings.h'), h_file.getvalue())

c_file = StringIO()
c_file.write("""/*
DO NOT EDIT!

This file is autogenerated by setup.py.  To edit its contents,
edit doc/docstrings.py

The weirdness here with strncpy is because some C compilers, notably
MSVC, do not support string literals greater than 256 characters.
*/

#include <string.h>
#include "docstrings.h"

""")
for key in keys:
    val = docstrings[key]
    c_file.write('char doc_%s[%d];\n' % (key, len(val)))

c_file.write("\nvoid fill_docstrings(void)\n{\n")
for key in keys:
    val = docstrings[key]
    # For portability across various compilers, we need to fill the
    # docstrings in 256-character chunks
    for i in range(0, len(val), 256):
        chunk = string_escape(val[i:i+256]).replace('"', '\\"')
        c_file.write('   strncpy(doc_%s + %d, "%s", %d);\n' % (
            key, i, chunk, min(len(val) - i, 256)))
    c_file.write("\n")
c_file.write("\n}\n\n")

write_if_different(join(srcroot, 'src', 'docstrings.c'), c_file.getvalue())

######################################################################
# PYWCS-SPECIFIC AND WRAPPER SOURCE FILES
PYWCS_VERSION = '1.11'
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
    'util.c',
    'wcslib_wrap.c',
    'wcslib_tabprm_wrap.c',
    'wcslib_units_wrap.c',
    'wcslib_wtbarr_wrap.c']
PYWCS_SOURCES = [join('src', x) for x in PYWCS_SOURCES]

######################################################################
# DISTUTILS SETUP

def get_distutils_option(option, commands):
    """ Returns the value of the given distutils option.

    Parameters
    ----------
    option : str
        The name of the option

    commands : list of str
        The list of commands on which this option is available

    Returns
    -------
    val : str or None
        the value of the given distutils option. If the option is not set,
        returns None.
    """
    # Pre-parse the Distutils command-line options and config files to
    # if the option is set.
    dist = Distribution()
    try:
        dist.parse_config_files()
        dist.parse_command_line()
    except DistutilsError:
        # Let distutils handle this itself
        return None
    except AttributeError:
        # This seems to get thrown for ./setup.py --help
        return None

    for cmd in commands:
        if cmd in dist.commands:
            break
    else:
        return None

    for cmd in commands:
        cmd_opts = dist.get_option_dict(cmd)
        if option in cmd_opts:
            return cmd_opts[option][1]
    else:
        return None


def adjust_compiler():
    """
    This function detects broken compilers and switches to another.  If
    the environment variable CC is explicitly set, or a compiler is
    specified on the commandline, no override is performed -- the purpose
    here is to only override a default compiler.

    The specific compilers with problems are:

        * The default compiler in XCode-4.2, llvm-gcc-4.2,
          segfaults when compiling wcslib.

    The set of broken compilers can be updated by changing the
    compiler_mapping variable.  It is a list of 2-tuples where the
    first in the pair is a regular expression matching the version
    of the broken compiler, and the second is the compiler to change
    to.
    """
    if 'CC' in os.environ:
        return

    if get_distutils_option(
        'compiler', ['build', 'build_ext', 'build_clib']) is not None:
        return

    from distutils import ccompiler
    import subprocess
    import re

    compiler_mapping = [
        ('i686-apple-darwin[0-9]*-llvm-gcc-4.2', 'clang')
        ]

    c = ccompiler.new_compiler()
    # The MSVC ccompiler class doesn't have a `compiler` member.
    if not hasattr(c, 'compiler'):
        return
    process = subprocess.Popen(
        c.compiler + ['--version'], stdout=subprocess.PIPE)
    output = process.communicate()[0].strip().decode('ascii')
    version = output.split()[0]
    for broken, fixed in compiler_mapping:
        if re.match(broken, version):
            os.environ['CC'] = fixed
            break


adjust_compiler()

libraries = []
define_macros = [('ECHO', None),
                 ('WCSTRIG_MACRO', None),
                 ('PYWCS_BUILD', None),
                 ('_GNU_SOURCE', None)]
undef_macros = []
extra_compile_args = []
extra_link_args = []
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
    # we are also writing these into wcsconfig.h, but that include
    # file is not used by wcslib
    define_macros.extend([
        ('YY_NO_UNISTD_H', None),
        ('_CRT_SECURE_NO_WARNINGS', None),
        ('_NO_OLDNAMES', None), # for mingw32
        ('NO_OLDNAMES', None), # for mingw64
        ('__STDC__', None) # for MSVC
        ])

if sys.platform.startswith('linux'):
    define_macros.append(('HAVE_SINCOS', None))

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
              extra_link_args=extra_link_args,
              libraries=libraries
              )
    ]

pkg = ["pywcs", "pywcs.tests"]

setupargs = {
    'version' :     VERSION,
    'description':  "Python wrappers to WCSLIB",
    'author' :      CONTACT,
    'author_email': EMAIL,
    'url' :         "http://projects.scipy.org/astropy/astrolib/wiki/WikiStart",
    'platforms' :   ["unix","windows"],
    'ext_modules' : PYWCS_EXTENSIONS,
    'data_files' : [
        ( 'pywcs/include', ['src/*.h']),
        ( 'pywcs/include/wcslib', [ WCSLIBC + '/*.h'] ),
        ( 'pywcs/tests/maps', ['lib/pywcs/tests/maps/*.hdr']),
        ( 'pywcs/tests/nrao', ['lib/pywcs/tests/nrao/*.hdr']),
        ( 'pywcs/tests/spectra', ['lib/pywcs/tests/spectra/*.hdr']),
        ( 'pywcs/tests/data', ['lib/pywcs/tests/data/*.hdr', 'lib/pywcs/tests/data/*.fits'])
        ],
    'package_dir' : { 'pywcs' : 'lib/pywcs', 'pywcs.tests' : 'lib/pywcs/tests'},
}


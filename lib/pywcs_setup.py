from __future__ import with_statement, division

import glob
import os
import shutil
import sys
import textwrap

from distutils.dist import Distribution 
from stsci.distutils.hooks import is_display_option


if sys.version_info[0] >= 3:
    def b(s):
        return s.encode('ascii')

    def string_escape(s):
        s = s.decode('ascii').encode('ascii', 'backslashreplace')
        s = s.replace(b('\n'), b('\\n'))
        return s.decode('ascii')
    from io import StringIO
else:
    def string_escape(s):
        return s.encode('string_escape')
    from cStringIO import StringIO


# BUILD should be 'debug', 'profile' or 'release'
# TODO: How often is this actually mucked with? Would it be worth adding a
# custom command that adds a command-line option for this?
BUILD = 'release'
OPENMP = False


def setup_hook(config):
    if is_display_option():
        return

    WCSLIB_VERSION = config['metadata']['stsci_wcs_version'].strip() 

    # Ensure that headers are copied to the correct location for installation
    # purposes
    if 'packages_root' in config['files']:
        root = config['files']['packages_root']
    else:
        root = ''
    include_dir = os.path.join(root, 'pywcs', 'include')
    wcslib_include_dir = os.path.join(include_dir, 'wcslib')
    for d in [include_dir, wcslib_include_dir]:
        if not os.path.exists(d):
            os.mkdir(d)
    for header in glob.glob(os.path.join('src', '*.h')):
        shutil.copy2(header, include_dir)
    for header in glob.glob(os.path.join('wcslib', 'C', '*.h')):
        shutil.copy2(header, wcslib_include_dir)

    # WCSLIB CONFIGURATION

    wcsconfig_h = StringIO()
    wcsconfig_h.write( wcsconfig_h_proto % 
            ( WCSLIB_VERSION, _determine_64_bit_int() )) 

    _write_if_different(os.path.join('src', 'wcsconfig.h'),
                        wcsconfig_h.getvalue())


    _generate_c_docstrings(config)

    # Some compilers don't work - work around that. 
    _adjust_compiler() 

    # Add/remove macros and compile args based on the build type
    libraries = []
    define_macros = []
    undef_macros = []
    extra_compile_args = []

    if BUILD.lower() == 'debug':
        define_macros.append(('DEBUG', None))
        undef_macros.append('NDEBUG')
        if not sys.platform.startswith('sun') and not sys.platform == 'win32':
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
        raise ValueError("BUILD should be one of 'debug', 'profile', or "
                         "'release'; got %s" % BUILD)

    if sys.platform == 'win32':
        # This duplicates a bunch of macros defined in wcsconfig_h_proto 
        # (below).  When compiling wcslib, not everything includes that 
        # include file, so we also define these macros on the command line. 

        define_macros.append(('YY_NO_UNISTD_H', None))
        define_macros.append(('_CRT_SECURE_NO_WARNINGS', None))
        define_macros.append(('_NO_OLDNAMES', None))
        define_macros.append(('NO_OLDNAMES', None))
        define_macros.append(('__STDC__', None))

    if sys.platform.startswith('linux'):
        define_macros.append(('HAVE_SINCOS', None))

    if not sys.platform.startswith('sun') and not sys.platform == 'win32':
        if OPENMP:
            extra_compile_args.append('-fopenmp')
            libraries.append('gomp')
        else:
            extra_compile_args.extend(['-Wno-unknown-pragmas'])

    for idx, m in enumerate(define_macros):
        if m[1] is not None:
            define_macros[idx] = '%s = %s' % m
        else:
            define_macros[idx] = m[0]

    ext_opts = [('libraries', libraries),
                ('define_macros', define_macros),
                ('undef_macros', undef_macros),
                ('extra_compile_args', extra_compile_args)]

    ext = config['extension=pywcs._pywcs']
    for opt, value in ext_opts:
        if opt in ext:
            ext[opt] += '\n' + '\n'.join(value)
        else:
            ext[opt] = '\n'.join(value)


def _generate_c_docstrings(config):
    # GENERATE DOCSTRINGS IN C
    docstrings = {}
    # We need to temporarily add lib/pywcs to the path for this to work
    if 'packages_root' in config['files']:
        root = config['files']['packages_root']
    else:
        root = ''
    pywcs_path = os.path.join(root, 'pywcs')
    if not pywcs_path in sys.path:
        clean_path = True
        sys.path.insert(0, pywcs_path)
    else:
        clean_path = False
    try:
        with open(os.path.join('doc', 'docstrings.py')) as f:
            exec(f.read(), docstrings)
    finally:
        if clean_path:
            sys.path.remove(pywcs_path)

    for key, val in list(docstrings.items()):
        if not (isinstance(key, str) and not key.startswith('__')):
            del docstrings[key]
    for key, val in docstrings.items():
        docstrings[key] = val.encode('utf8').lstrip()

    docstrings = sorted(docstrings.items())

    docstrings_h = StringIO()
    docstrings_h.write(textwrap.dedent("""
        /*
        DO NOT EDIT!

        This file is autogenerated by setup.p.  To edit its contents, edit
        doc/docstrings.py
        */

        #ifndef __DOCSTRINGS_H__
        #define __DOCSTRINGS_H__

        void fill_docstrings(void);

    """))

    for key, val in docstrings:
        docstrings_h.write('extern char doc_%s[%d];\n' % (key, len(val)))
    docstrings_h.write('\n#endif\n\n')

    _write_if_different(os.path.join('src', 'docstrings.h'),
                        docstrings_h.getvalue())

    docstrings_c = StringIO()
    docstrings_c.write(textwrap.dedent("""
        /*
        DO NOT EDIT!

        This file is autogenerated by setup.py.  To edit its contents,
        edit doc/docstrings.py

        The weirdness here with strncpy is because some C compilers, notably
        MSVC, do not support string literals greater than 256 characters.
        */

        #include <string.h>
        #include "docstrings.h"
    """))

    for key, val in docstrings:
        docstrings_c.write('char doc_%s[%d];\n' % (key, len(val)))

    docstrings_c.write('\nvoid fill_docstrings(void)\n{\n')
    for key, val in docstrings:
        # For portability across various compilers, we need to fill the
        # docstrings in 256-character chunks
        for idx in range(0, len(val), 256):
            chunk = string_escape(val[idx:idx + 256]).replace('"', '\\"')
            docstrings_c.write('   strncpy(doc_%s + %d, "%s", %d);\n'
                               % (key, idx, chunk, min(len(val) - idx, 256)))
        docstrings_c.write('\n')
    docstrings_c.write('\n}\n\n')

    _write_if_different(os.path.join('src', 'docstrings.c'),
                        docstrings_c.getvalue())


# The only configuration parameter needed at compile-time is how to 
# specify a 64-bit signed integer.  Python's ctypes module can get us 
# that information, but it is only available in Python 2.5 or later. 
# If we can't be absolutely certain, we default to "long long int", 
# which is correct on most platforms (x86, x86_64).  If we find 
# platforms where this heuristic doesn't work, we may need to hardcode 
# for them. 

def _determine_64_bit_int():
    try:
        try:
            import ctypes
        except ImportError:
            raise ValueError

        if ctypes.sizeof(ctypes.c_longlong) == 8:
            return 'long long int'
        elif ctypes.sizeof(ctypes.c_long) == 8:
            return 'long int'
        elif ctypes.sizeof(ctypes.c_int) == 8:
            return 'int'
        else:
            raise ValueError
    except ValueError:
        return 'long long int'


def _write_if_different(filename, data):
    data = data.encode('ascii')

    if os.path.exists(filename):
        fd = open(filename, 'rb')
        original_data = fd.read()
        fd.close()
    else:
        original_data = None

    if original_data != data:
        fd = open(filename, 'wb')
        fd.write(data)
        fd.close()

#

def _adjust_compiler():
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


# The prototype for the standard header file

wcsconfig_h_proto = """
/* WCSLIB library version number. */
#define WCSLIB_VERSION %s

/* 64-bit integer data type. */
#define WCSLIB_INT64 %s

/* Windows needs some other defines */
/* I think the only one we need is _WIN32, but the others don't hurt - Mark S. 2012-02-14 */
#if defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__) || defined (__MINGW64__)

/* 
* There is a lexical analyzer that was generated by flex.  Tell it
* that we don't have unistd.h
*/
#ifndef YY_NO_UNISTD_H
#define YY_NO_UNISTD_H
#endif

/*
* Visual C++ will warn you about certain functions that have a
* high risk of buffer overflows, such as strcpy().  This
* tells the compiler to shut up about it.
*/
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

/*
* Suppress definition of wcsset() by system include files.
* I don't know what's up with mingw.
* In Visual C, the __STDC__ excludes non-standard things that
* would otherwise be defined in the header files.
*/

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
"""


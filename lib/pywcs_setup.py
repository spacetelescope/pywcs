from __future__ import with_statement, division

import glob
import os
import shutil
import sys
import textwrap

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
WCSLIB_VERSION = '4.7' # TODO: Could this be determined automagically?


def setup_hook(config):
    if is_display_option():
        return

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
    # The only configuration parameter needed at compile-time is how to
    # specify a 64-bit signed integer.  Python's ctypes module can get us
    # that information, but it is only available in Python 2.5 or later.
    # If we can't be absolutely certain, we default to "long long int",
    # which is correct on most platforms (x86, x86_64).  If we find
    # platforms where this heuristic doesn't work, we may need to hardcode
    # for them.

    wcsconfig_h = StringIO()
    wcsconfig_h.write(textwrap.dedent("""
        /* WCSLIB library version number. */
        #define WCSLIB_VERSION %s

        /* 64-bit integer data type. */
        #define WCSLIB_INT64 %s
     """ % (WCSLIB_VERSION, _determine_64_bit_int())))

    _write_if_different(os.path.join('src', 'wcsconfig.h'),
                        wcsconfig_h.getvalue())


    _generate_c_docstrings(config)


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
        define_macros.append(('YY_NO_UNISTD_H', None))
        define_macros.append(('_CRT_SECURE_NO_WARNINGS', None))
        define_macros.append(('_NO_OLDNAMES', None)) # for mingw32
        define_macros.append(('NO_OLDNAMES', None)) # for mingw64

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

    for key, val in docstrings.items():
        if not (isinstance(key, basestring) and not key.startswith('__')):
            del docstrings[key]
    for key, val in docstrings.items():
        docstrings[key] = val.encode('utf8').lstrip()

    docstrings = sorted(docstrings.iteritems())

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

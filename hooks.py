import os
from distutils.dist import Distribution

def pywcs_hook(config) :

    WCSVERSION = config['metadata']['stsci_wcs_version'].strip()

    ## write wcsconfig.h
    write_if_different( 'src/wcsconfig.h', wcsconfig_h_proto % ( WCSVERSION, determine_64_bit_int()) )

    ## circumvent buggy compilers
    adjust_compiler()

    return
    

######################################################################

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

######################################################################u

wcsconfig_h_proto = """
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
"""

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


######################################################################
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



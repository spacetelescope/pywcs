def _fix(content, indent=0):
    lines = content.split('\n')
    indent = '\n' + ' ' * indent
    return indent.join(lines)

def ONE_OR_TWO_ARGS(out_type, indent=0):
    return _fix("""Either one or two arguments may be provided.

    - one argument: An Nx2 array of I{x}- and I{y}-coordinates.

    - two arguments: Two one-dimensional arrays of I{x} and I{y}
          coordinates.

@return: Returns the %s coordinates.  If the input was a
     single array, a single array is returned, otherwise a tuple of
     arrays is returned.""" % out_type, indent)

def FITS_EQUIVALENT(method_name, indent=0):
    return _fix("""B{The pixel coordinates are 0-based (like array indices in C and
Python).  If your pixel coordinates are 1-based (like array
indices in Fortran), use L{%s_fits} instead.}""" % method_name, indent)

def NON_FITS_EQUIVALENT(method_name, indent=0):
    return _fix("""Identical to L{%s}, except pixel coordinates are 1-based (like array
indices in Fortran), instead of 0-based (like array indices C and
Python).""" % method_name, indent)


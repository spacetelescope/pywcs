#!/usr/bin/env python
from __future__ import division # confidence high

try:
    import stsci.tools.stsci_distutils_hack as H
except ImportError:
    import stsci_distutils_hack as H
H.run()

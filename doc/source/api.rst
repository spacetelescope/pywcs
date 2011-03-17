.. include:: references.rst

API documentation
=================

:mod:`pywcs`
------------

.. automodule:: pywcs.pywcs

Classes
-------

.. toctree::
   :maxdepth: 2

   api_wcs.rst
   api_wcsprm.rst
   api_distortion.rst
   api_sip.rst
   api_units.rst
   relax.rst

Testing pywcs
=============

The unit tests are written for use with `nose
<http://code.google.com/p/python-nose/>`.  To run them, install pywcs
and then at the commandline::

   nosetests pywcs.tests


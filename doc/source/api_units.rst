.. include:: references.rst

`UnitConverter`
```````````````

.. autoclass:: pywcs.UnitConverter
   :members:
   :inherited-members:
   :undoc-members:

.. _unit-aliases:

Unit aliases
------------

When converting non-standard units to standard ones, a case-sensitive
match is required for the aliases listed below, in particular the only
recognized aliases with metric prefixes are ``"KM"``, ``"KHZ"``,
``"MHZ"``, and ``"GHZ"``.  Potentially unsafe translations of ``"D"``,
``"H"``, and ``"S"``, shown in parentheses, are optional, using the
*translate_units* parameter.

========== =============================================================
Unit       Recognized aliases
========== =============================================================
Angstrom   angstrom
arcmin     arcmins, ARCMIN, ARCMINS
arcsec     arcsecs, ARCSEC, ARCSECS
beam       BEAM
byte       Byte
d          day, days, (D), DAY, DAYS
deg        degree, degrees, DEG, DEGREE, DEGREES
GHz        GHZ
h          hr, (H), HR
Hz         hz, HZ
kHz        KHZ
Jy         JY
K          kelvin, kelvins, Kelvin, Kelvins, KELVIN, KELVINS
km         KM
m          metre, meter, metres, meters, M, METRE, METER, METRES, METERS
min        MIN
MHz        MHZ
Ohm        ohm
Pa         pascal, pascals, Pascal, Pascals, PASCAL, PASCALS
pixel      pixels, PIXEL, PIXELS
rad        radian, radians, RAD, RADIAN, RADIANS
s          sec, second, seconds, (S), SEC, SECOND, SECONDS
V          volt, volts, Volt, Volts, VOLT, VOLTS
yr         year, years, YR, YEAR, YEARS
========== =============================================================

The aliases ``"angstrom"``, ``"ohm"``, and ``"Byte"`` for (Angstrom,
Ohm, and byte) are recognized by pywcs/wcslib itself as an unofficial
extension of the standard, but they are converted to the standard form
here.

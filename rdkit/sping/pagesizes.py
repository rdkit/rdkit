"""pagesizes.py

This module defines a few common page sizes in PIDDLE units (1/72 inch).
"""

from rdkit.sping.pid import cm, inch

_W, _H = (21 * cm, 29.7 * cm)

A6 = (_W * .5, _H * .5)
A5 = (_H * .5, _W)
A4 = (_W, _H)
A3 = (_H * 2, _W)
A2 = (_W * 2, _H * 2)
A1 = (_H * 4, _W * 2)
A0 = (_W * 4, _H * 4)

letter = (8.5 * inch, 11 * inch)
legal = (8.5 * inch, 17 * inch)
elevenSeventeen = (11 * inch, 17 * inch)

_BW, _BH = (25 * cm, 35.3 * cm)
B6 = (_BW * .5, _BH * .5)
B5 = (_BH * .5, _BW)
B4 = (_BW, _BH)
B3 = (_BH * 2, _BW)
B2 = (_BW * 2, _BH * 2)
B1 = (_BH * 4, _BW * 2)
B0 = (_BW * 4, _BH * 4)

#!/usr/bin/env python

lines = []
lines.append("This is a test of the <i>stringformat</i> module.")
lines.append("This module allows you to insert html-like formating tags")
lines.append("into strings. Built on top of sping and xmllib, stringformat.py")

lines.append(
  "allows use of <b>boldface</b>, <i>italic</i>, <u>underline</u>, <super>super</super>script,")
lines.append("<sub>sub</sub>script and greek letter symbols like &alpha;, &omega;, and &phi;, as")
lines.append("specified in MathML")

print("""
This is an example using formatted strings.  It renders its output to
a portable document file called "formatted-strings.pdf" and a
postscript file called "formatted-strings.ps" for viewing with Acrobat
Reader or for printing to a postscript printer.
""")

import sping.stringformat
from sping.PDF import PDFCanvas
from sping.PS import PSCanvas

# Do PDF first
canvas = PDFCanvas(size=(350, 200), name="formatted-strings.pdf")

y = 20
for line in lines:
  sping.stringformat.drawString(canvas, line, 10, y)
  y = y + 20

canvas.flush()
canvas.save()

# Now do postscript
canvas = PSCanvas(size=(350, 200), name="formatted-strings.ps")

y = 20
for line in lines:
  sping.stringformat.drawString(canvas, line, 10, y)
  y = y + 20

canvas.flush()
canvas.save()

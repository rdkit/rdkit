# $Id$
from __future__ import print_function
import pidtest


def testLatin1Chars(can):
  curx, cury = 10, 20
  can.drawString("hola M�laga amigos ni�os", curx, cury)
  cury = cury + 20
  can.drawString("Como estan?", curx, cury)

  can.flush()
  can.clear()  # get our next page ???
  curx, cury = 10, 20
  can.drawString("hola M�laga amigos ni�os: Page 2", curx, cury)

  str = "sometext��� with �c�nts"
  print(len("���"))
  pidtest.CenterAndBox(can, str, y=150)

  can.flush()

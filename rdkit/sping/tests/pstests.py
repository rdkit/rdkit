# -*- coding: utf-8 -*-

import pidtest


def testLatin1Chars(can):
  curx, cury = 10, 20
  can.drawString(u"años luz detrás", curx, cury)
  cury = cury + 20
  can.drawString("Como estan?", curx, cury)

  can.flush()
  can.clear()  # get our next page ???
  curx, cury = 10, 20
  can.drawString(u"años luz detrás: page 2", curx, cury)

  str = u"some text with ácceñts"
  pidtest.CenterAndBox(can, str, y=150)

  can.flush()

# $Id: pstests.py 3194 2004-02-24 23:58:17Z glandrum $
import pidtest

def testLatin1Chars(can):
    curx,cury = 10,20
    can.drawString("hola Málaga amigos niños", curx,cury)
    cury = cury + 20
    can.drawString("Como estan?", curx, cury)

    can.flush()
    can.clear()  # get our next page ???
    curx,cury = 10,20
    can.drawString("hola Málaga amigos niños: Page 2", curx,cury)

    str = "sometextÄËÖ with ácënts"
    print len("ÄËÖ")
    pidtest.CenterAndBox(can, str, y=150)

    can.flush()



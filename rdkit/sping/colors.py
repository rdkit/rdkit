# -------------------------------------------------------------------------
# Color
# -------------------------------------------------------------------------


class Color:
  """This class is used to represent color.  Components red, green, blue 
          are in the range 0 (dark) to 1 (full intensity)."""

  def __init__(self, red=0, green=0, blue=0):
    "Initialize with red, green, blue in range [0-1]."
    _float = float
    d = self.__dict__
    d["red"] = _float(red)
    d["green"] = _float(green)
    d["blue"] = _float(blue)

  def __setattr__(self, name, value):
    raise TypeError("piddle.Color has read-only attributes")

  def __mul__(self, x):
    return Color(self.red * x, self.green * x, self.blue * x)

  def __rmul__(self, x):
    return Color(self.red * x, self.green * x, self.blue * x)

  def __truediv__(self, x):
    return Color(self.red / x, self.green / x, self.blue / x)

  def __div__(self, x):
    return Color(self.red / x, self.green / x, self.blue / x)

  def __rdiv__(self, x):
    return Color(self.red / x, self.green / x, self.blue / x)

  def __add__(self, x):
    return Color(self.red + x.red, self.green + x.green, self.blue + x.blue)

  def __sub__(self, x):
    return Color(self.red - x.red, self.green - x.green, self.blue - x.blue)

  def __repr__(self):
    return "Color(%1.2f,%1.2f,%1.2f)" % (self.red, self.green, self.blue)

  def __hash__(self):
    return hash((self.red, self.green, self.blue))

  def __cmp__(self, other):
    try:
      dsum = 4 * self.red - 4 * other.red + 2 * self.green - 2 * other.green + self.blue - other.blue
    except Exception:
      return -1
    if dsum > 0:
      return 1
    if dsum < 0:
      return -1
    return 0

  def toHexRGB(self):
    "Convert the color back to an integer suitable for the "
    "0xRRGGBB hex representation"
    r = int(0xFF * self.red)
    g = int(0xFF * self.green)
    b = int(0xFF * self.blue)
    # print "r= %d, g=%d, b = %d" % (r,b,g)
    return (r << 16) + (g << 8) + b

  def toHexStr(self):
    return "0x%.6x" % self.toHexRGB()


def HexColor(val):
  """This class converts a hex string, or an actual integer number,
          into the corresponding color.  E.g., in "AABBCC" or 0xAABBCC,
          AA is the red, BB is the green, and CC is the blue (00-FF)."""
  if isinstance(val, str):
    val = int(val, 16)
  factor = 1.0 / 255
  return Color(factor * ((val >> 16) & 0xFF), factor * ((val >> 8) & 0xFF), factor * (val & 0xFF))


# color constants -- mostly from HTML standard
aliceblue = HexColor(0xF0F8FF)
antiquewhite = HexColor(0xFAEBD7)
aqua = HexColor(0x00FFFF)
aquamarine = HexColor(0x7FFFD4)
azure = HexColor(0xF0FFFF)
beige = HexColor(0xF5F5DC)
bisque = HexColor(0xFFE4C4)
black = HexColor(0x000000)
blanchedalmond = HexColor(0xFFEBCD)
blue = HexColor(0x0000FF)
blueviolet = HexColor(0x8A2BE2)
brown = HexColor(0xA52A2A)
burlywood = HexColor(0xDEB887)
cadetblue = HexColor(0x5F9EA0)
chartreuse = HexColor(0x7FFF00)
chocolate = HexColor(0xD2691E)
coral = HexColor(0xFF7F50)
cornflower = HexColor(0x6495ED)
cornsilk = HexColor(0xFFF8DC)
crimson = HexColor(0xDC143C)
cyan = HexColor(0x00FFFF)
darkblue = HexColor(0x00008B)
darkcyan = HexColor(0x008B8B)
darkgoldenrod = HexColor(0xB8860B)
darkgray = HexColor(0xA9A9A9)
darkgreen = HexColor(0x006400)
darkkhaki = HexColor(0xBDB76B)
darkmagenta = HexColor(0x8B008B)
darkolivegreen = HexColor(0x556B2F)
darkorange = HexColor(0xFF8C00)
darkorchid = HexColor(0x9932CC)
darkred = HexColor(0x8B0000)
darksalmon = HexColor(0xE9967A)
darkseagreen = HexColor(0x8FBC8B)
darkslateblue = HexColor(0x483D8B)
darkslategray = HexColor(0x2F4F4F)
darkturquoise = HexColor(0x00CED1)
darkviolet = HexColor(0x9400D3)
deeppink = HexColor(0xFF1493)
deepskyblue = HexColor(0x00BFFF)
dimgray = HexColor(0x696969)
dodgerblue = HexColor(0x1E90FF)
firebrick = HexColor(0xB22222)
floralwhite = HexColor(0xFFFAF0)
forestgreen = HexColor(0x228B22)
fuchsia = HexColor(0xFF00FF)
gainsboro = HexColor(0xDCDCDC)
ghostwhite = HexColor(0xF8F8FF)
gold = HexColor(0xFFD700)
goldenrod = HexColor(0xDAA520)
gray = HexColor(0x808080)
grey = gray
green = HexColor(0x008000)
greenyellow = HexColor(0xADFF2F)
honeydew = HexColor(0xF0FFF0)
hotpink = HexColor(0xFF69B4)
indianred = HexColor(0xCD5C5C)
indigo = HexColor(0x4B0082)
ivory = HexColor(0xFFFFF0)
khaki = HexColor(0xF0E68C)
lavender = HexColor(0xE6E6FA)
lavenderblush = HexColor(0xFFF0F5)
lawngreen = HexColor(0x7CFC00)
lemonchiffon = HexColor(0xFFFACD)
lightblue = HexColor(0xADD8E6)
lightcoral = HexColor(0xF08080)
lightcyan = HexColor(0xE0FFFF)
lightgoldenrodyellow = HexColor(0xFAFAD2)
lightgreen = HexColor(0x90EE90)
lightgrey = HexColor(0xD3D3D3)
lightpink = HexColor(0xFFB6C1)
lightsalmon = HexColor(0xFFA07A)
lightseagreen = HexColor(0x20B2AA)
lightskyblue = HexColor(0x87CEFA)
lightslategray = HexColor(0x778899)
lightsteelblue = HexColor(0xB0C4DE)
lightyellow = HexColor(0xFFFFE0)
lime = HexColor(0x00FF00)
limegreen = HexColor(0x32CD32)
linen = HexColor(0xFAF0E6)
magenta = HexColor(0xFF00FF)
maroon = HexColor(0x800000)
mediumaquamarine = HexColor(0x66CDAA)
mediumblue = HexColor(0x0000CD)
mediumorchid = HexColor(0xBA55D3)
mediumpurple = HexColor(0x9370DB)
mediumseagreen = HexColor(0x3CB371)
mediumslateblue = HexColor(0x7B68EE)
mediumspringgreen = HexColor(0x00FA9A)
mediumturquoise = HexColor(0x48D1CC)
mediumvioletred = HexColor(0xC71585)
midnightblue = HexColor(0x191970)
mintcream = HexColor(0xF5FFFA)
mistyrose = HexColor(0xFFE4E1)
moccasin = HexColor(0xFFE4B5)
navajowhite = HexColor(0xFFDEAD)
navy = HexColor(0x000080)
oldlace = HexColor(0xFDF5E6)
olive = HexColor(0x808000)
olivedrab = HexColor(0x6B8E23)
orange = HexColor(0xFFA500)
orangered = HexColor(0xFF4500)
orchid = HexColor(0xDA70D6)
palegoldenrod = HexColor(0xEEE8AA)
palegreen = HexColor(0x98FB98)
paleturquoise = HexColor(0xAFEEEE)
palevioletred = HexColor(0xDB7093)
papayawhip = HexColor(0xFFEFD5)
peachpuff = HexColor(0xFFDAB9)
peru = HexColor(0xCD853F)
pink = HexColor(0xFFC0CB)
plum = HexColor(0xDDA0DD)
powderblue = HexColor(0xB0E0E6)
purple = HexColor(0x800080)
red = HexColor(0xFF0000)
rosybrown = HexColor(0xBC8F8F)
royalblue = HexColor(0x4169E1)
saddlebrown = HexColor(0x8B4513)
salmon = HexColor(0xFA8072)
sandybrown = HexColor(0xF4A460)
seagreen = HexColor(0x2E8B57)
seashell = HexColor(0xFFF5EE)
sienna = HexColor(0xA0522D)
silver = HexColor(0xC0C0C0)
skyblue = HexColor(0x87CEEB)
slateblue = HexColor(0x6A5ACD)
slategray = HexColor(0x708090)
snow = HexColor(0xFFFAFA)
springgreen = HexColor(0x00FF7F)
steelblue = HexColor(0x4682B4)
tan = HexColor(0xD2B48C)
teal = HexColor(0x008080)
thistle = HexColor(0xD8BFD8)
tomato = HexColor(0xFF6347)
turquoise = HexColor(0x40E0D0)
violet = HexColor(0xEE82EE)
wheat = HexColor(0xF5DEB3)
white = HexColor(0xFFFFFF)
whitesmoke = HexColor(0xF5F5F5)
yellow = HexColor(0xFFFF00)
yellowgreen = HexColor(0x9ACD32)

# special case -- indicates no drawing should be done
transparent = Color(-1, -1, -1)

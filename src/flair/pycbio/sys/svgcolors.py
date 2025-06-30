# Copyright 2006-2025 Mark Diekhans

# table lifted from
#  https://pypi.org/project/webcolors/

from flair.pycbio import PycbioException
from flair.pycbio.sys.color import Color

def _mkcolor(rgb8):
    return Color.fromPackRgb8(rgb8)


class SvgColors:
    """Definitions of SVG/CCS3 colors by name
    good table here: http://www.december.com/html/spec/colorsvg.html.
    Some have duplicate names.
    """

    aliceblue = _mkcolor(0xf0f8ff)
    antiquewhite = _mkcolor(0xfaebd7)
    aqua = _mkcolor(0x00ffff)
    aquamarine = _mkcolor(0x7fffd4)
    azure = _mkcolor(0xf0ffff)
    beige = _mkcolor(0xf5f5dc)
    bisque = _mkcolor(0xffe4c4)
    black = _mkcolor(0x000000)
    blanchedalmond = _mkcolor(0xffebcd)
    blue = _mkcolor(0x0000ff)
    blueviolet = _mkcolor(0x8a2be2)
    brown = _mkcolor(0xa52a2a)
    burlywood = _mkcolor(0xdeb887)
    cadetblue = _mkcolor(0x5f9ea0)
    chartreuse = _mkcolor(0x7fff00)
    chocolate = _mkcolor(0xd2691e)
    coral = _mkcolor(0xff7f50)
    cornflowerblue = _mkcolor(0x6495ed)
    cornsilk = _mkcolor(0xfff8dc)
    crimson = _mkcolor(0xdc143c)
    cyan = _mkcolor(0x00ffff)
    darkblue = _mkcolor(0x00008b)
    darkcyan = _mkcolor(0x008b8b)
    darkgoldenrod = _mkcolor(0xb8860b)
    darkgray = _mkcolor(0xa9a9a9)
    darkgrey = _mkcolor(0xa9a9a9)
    darkgreen = _mkcolor(0x006400)
    darkkhaki = _mkcolor(0xbdb76b)
    darkmagenta = _mkcolor(0x8b008b)
    darkolivegreen = _mkcolor(0x556b2f)
    darkorange = _mkcolor(0xff8c00)
    darkorchid = _mkcolor(0x9932cc)
    darkred = _mkcolor(0x8b0000)
    darksalmon = _mkcolor(0xe9967a)
    darkseagreen = _mkcolor(0x8fbc8f)
    darkslateblue = _mkcolor(0x483d8b)
    darkslategray = _mkcolor(0x2f4f4f)
    darkslategrey = _mkcolor(0x2f4f4f)
    darkturquoise = _mkcolor(0x00ced1)
    darkviolet = _mkcolor(0x9400d3)
    deeppink = _mkcolor(0xff1493)
    deepskyblue = _mkcolor(0x00bfff)
    dimgray = _mkcolor(0x696969)
    dimgrey = _mkcolor(0x696969)
    dodgerblue = _mkcolor(0x1e90ff)
    firebrick = _mkcolor(0xb22222)
    floralwhite = _mkcolor(0xfffaf0)
    forestgreen = _mkcolor(0x228b22)
    fuchsia = _mkcolor(0xff00ff)
    gainsboro = _mkcolor(0xdcdcdc)
    ghostwhite = _mkcolor(0xf8f8ff)
    gold = _mkcolor(0xffd700)
    goldenrod = _mkcolor(0xdaa520)
    gray = _mkcolor(0x808080)
    grey = gray
    green = _mkcolor(0x008000)
    greenyellow = _mkcolor(0xadff2f)
    honeydew = _mkcolor(0xf0fff0)
    hotpink = _mkcolor(0xff69b4)
    indianred = _mkcolor(0xcd5c5c)
    indigo = _mkcolor(0x4b0082)
    ivory = _mkcolor(0xfffff0)
    khaki = _mkcolor(0xf0e68c)
    lavender = _mkcolor(0xe6e6fa)
    lavenderblush = _mkcolor(0xfff0f5)
    lawngreen = _mkcolor(0x7cfc00)
    lemonchiffon = _mkcolor(0xfffacd)
    lightblue = _mkcolor(0xadd8e6)
    lightcoral = _mkcolor(0xf08080)
    lightcyan = _mkcolor(0xe0ffff)
    lightgoldenrodyellow = _mkcolor(0xfafad2)
    lightgrey = _mkcolor(0xd3d3d3)
    lightgreen = _mkcolor(0x90ee90)
    lightpink = _mkcolor(0xffb6c1)
    lightsalmon = _mkcolor(0xffa07a)
    lightseagreen = _mkcolor(0x20b2aa)
    lightskyblue = _mkcolor(0x87cefa)
    lightslategrey = _mkcolor(0x778899)
    lightsteelblue = _mkcolor(0xb0c4de)
    lightyellow = _mkcolor(0xffffe0)
    lime = _mkcolor(0x00ff00)
    limegreen = _mkcolor(0x32cd32)
    linen = _mkcolor(0xfaf0e6)
    magenta = _mkcolor(0xff00ff)
    maroon = _mkcolor(0x800000)
    mediumaquamarine = _mkcolor(0x66cdaa)
    mediumblue = _mkcolor(0x0000cd)
    mediumorchid = _mkcolor(0xba55d3)
    mediumpurple = _mkcolor(0x9370d8)
    mediumseagreen = _mkcolor(0x3cb371)
    mediumslateblue = _mkcolor(0x7b68ee)
    mediumspringgreen = _mkcolor(0x00fa9a)
    mediumturquoise = _mkcolor(0x48d1cc)
    mediumvioletred = _mkcolor(0xc71585)
    midnightblue = _mkcolor(0x191970)
    mintcream = _mkcolor(0xf5fffa)
    mistyrose = _mkcolor(0xffe4e1)
    moccasin = _mkcolor(0xffe4b5)
    navajowhite = _mkcolor(0xffdead)
    navy = _mkcolor(0x000080)
    oldlace = _mkcolor(0xfdf5e6)
    olive = _mkcolor(0x808000)
    olivedrab = _mkcolor(0x6b8e23)
    orange = _mkcolor(0xffa500)
    orangered = _mkcolor(0xff4500)
    orchid = _mkcolor(0xda70d6)
    palegoldenrod = _mkcolor(0xeee8aa)
    palegreen = _mkcolor(0x98fb98)
    paleturquoise = _mkcolor(0xafeeee)
    palevioletred = _mkcolor(0xd87093)
    papayawhip = _mkcolor(0xffefd5)
    peachpuff = _mkcolor(0xffdab9)
    peru = _mkcolor(0xcd853f)
    pink = _mkcolor(0xffc0cb)
    plum = _mkcolor(0xdda0dd)
    powderblue = _mkcolor(0xb0e0e6)
    purple = _mkcolor(0x800080)
    red = _mkcolor(0xff0000)
    rosybrown = _mkcolor(0xbc8f8f)
    royalblue = _mkcolor(0x4169e1)
    saddlebrown = _mkcolor(0x8b4513)
    salmon = _mkcolor(0xfa8072)
    sandybrown = _mkcolor(0xf4a460)
    seagreen = _mkcolor(0x2e8b57)
    seashell = _mkcolor(0xfff5ee)
    sienna = _mkcolor(0xa0522d)
    silver = _mkcolor(0xc0c0c0)
    skyblue = _mkcolor(0x87ceeb)
    slateblue = _mkcolor(0x6a5acd)
    slategray = _mkcolor(0x708090)
    slategrey = _mkcolor(0x708090)
    snow = _mkcolor(0xfffafa)
    springgreen = _mkcolor(0x00ff7f)
    steelblue = _mkcolor(0x4682b4)
    tan = _mkcolor(0xd2b48c)
    teal = _mkcolor(0x008080)
    thistle = _mkcolor(0xd8bfd8)
    tomato = _mkcolor(0xff6347)
    turquoise = _mkcolor(0x40e0d0)
    violet = _mkcolor(0xee82ee)
    wheat = _mkcolor(0xf5deb3)
    white = _mkcolor(0xffffff)
    whitesmoke = _mkcolor(0xf5f5f5)
    yellow = _mkcolor(0xffff00)
    yellowgreen = _mkcolor(0x9acd32)

    # Parallel arrays with names and colors, sorted by name.
    # Built in a lazy manner
    _names = None
    _colors = None

    # table of id of color objects to name
    _svgcolorIdNameMap = {}

    @classmethod
    def _initColors(cls):
        names = []
        colors = []
        for n in sorted(vars(cls).keys()):
            c = getattr(cls, n)
            if isinstance(c, Color):
                names.append(n)
                colors.append(c)
                cls._svgcolorIdNameMap[id(c)] = n
        cls._names = tuple(names)
        cls._colors = tuple(colors)

    @classmethod
    def lookup(cls, name):
        "convert case-insensitive name to color"
        try:
            return getattr(cls, name.lower())
        except AttributeError:
            raise PycbioException(f"unknown SVG color '{name}'")

    @classmethod
    def getName(cls, color):
        """Get the name for a color. The color must be one of the members
        of this class"""
        if cls._names is None:
            cls._initColors()
        try:
            return cls._svgcolorIdNameMap[id(color)]
        except KeyError:
            if not isinstance(color, Color):
                raise PycbioException(f"object is not a Color '{type(color)}'")
            else:
                raise PycbioException(f"color is not an SVG color object '{color.toRgbStr()}'")

    @classmethod
    def getNames(cls):
        "get all color names"
        if cls._names is None:
            cls._initColors()
        return cls._names

    @classmethod
    def _findClosestIdx(cls, color):
        if cls._names is None:
            cls._initColors()
        minIdx = 0
        minDist = color.distance(cls._colors[minIdx])
        for i in range(1, len(cls._colors)):
            d = color.distance(cls._colors[i])
            if d < minDist:
                minIdx = i
                minDist = d
        return minIdx

    @classmethod
    def getClosestName(cls, color):
        "get the closet SVG name"
        i = cls._findClosestIdx(color)
        return cls._names[i]

    @classmethod
    def getClosestColor(cls, color):
        "get the closet SVG name"
        i = cls._findClosestIdx(color)
        return cls._colors[i]

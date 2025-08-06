# Copyright 2006-2025 Mark Diekhans
import re
import colorsys
import math

# FIXME: do we really need to store both RGB and HSV?, maybe cache HSV converson


def _int8ToReal(v):
    assert 0 <= v <= 255
    return v / 255.0

def _realToInt8(v):
    return int(round(v * 255))

class Color:
    """Immutable color object with conversion to/from different formats.
    Don't construct directly use, factory (from*) static methods.

    :ivar red: red channel, in the range 0.0..1.0
    :ivar green: green channel, in the range 0.0..1.0
    :ivar blue: blue channel, in the range 0.0..1.0
    :ivar alpha: it not None, alpha value in the range 0.0..1.0
    """
    __slots__ = ("red", "green", "blue", "alpha", "_hsv")

    def __init__(self, red, green, blue, alpha=None, *, hsv=None):
        object.__setattr__(self, 'red', red)
        object.__setattr__(self, 'green', green)
        object.__setattr__(self, 'blue', blue)
        object.__setattr__(self, 'alpha', alpha)
        object.__setattr__(self, '_hsv', None)

    def __getstate__(self):
        return (self.red, self.green, self.blue, self.alpha)

    def __setstate__(self, state):
        object.__setattr__(self, 'red', state[0])
        object.__setattr__(self, 'green', state[1])
        object.__setattr__(self, 'blue', state[2])
        object.__setattr__(self, 'alpha', state[3])
        object.__setattr__(self, '_hsv', None)

    def __str__(self):
        if self.alpha is None:
            return str((self.red, self.green, self.blue))
        else:
            return str((self.red, self.green, self.blue, self.alpha))

    def __repr__(self):
        return f"Color({self.red}, {self.green}, {self.blue}, {self.alpha}"

    def __setattr__(self, key, value):
        raise AttributeError(f"{self.__class__.__name__} is immutable")

    def __hash__(self):
        return hash((self.red, self.green, self.blue, self.alpha))

    def __eq__(self, other):
        if isinstance(other, Color):
            return (self.red, self.green, self.blue, self.alpha) == (other.red, other.green, other.blue, other.alpha)
        return False

    def _createHsvCache(self):
        "create HSV cache;  idempotent so thread-safe"
        object.__setattr__(self, '_hsv',
                           colorsys.rgb_to_hsv(self.red, self.green, self.blue))

    @property
    def hue(self):
        "get hue component, in the range 0.0..1.0"
        if self._hsv is None:
            self._createHsvCache()
        return self._hsv[0]

    @property
    def saturation(self):
        "get saturation component, in the range 0.0..1.0"
        if self._hsv is None:
            self._createHsvCache()
        return self._hsv[1]

    @property
    def value(self):
        "get value component, in the range 0.0..1.0"
        if self._hsv is None:
            self._createHsvCache()
        return self._hsv[2]

    @property
    def alphaDflt(self):
        "alpha, defaulted to 1.0 if None"
        return self.alpha if self.alpha is not None else 1.0

    @property
    def red8(self):
        "red channel as an 8-bit int"
        return _realToInt8(self.red)

    @property
    def green8(self):
        "green channel as an 8-bit int"
        return _realToInt8(self.green)

    @property
    def blue8(self):
        "get blue channel as an 8-bit int"
        return _realToInt8(self.blue)

    @property
    def alpha8(self):
        "get alpha as an 8-bit int, or None if no alpha"
        return _realToInt8(self.alpha) if self.alpha is not None else None

    @property
    def hue8(self):
        "hue component as an integer angle"
        return int(round(360 * self.hue))

    @property
    def saturation8(self):
        "saturation component as an integer percent"
        return int(round(100 * self.saturation))

    @property
    def value8(self):
        "value component as an integer percent"
        return int(round(100 * self.value))

    @property
    def rgb(self):
        "RGB as tuple of real numbers"
        return (self.red, self.green, self.blue)

    @property
    def rgba(self):
        "RGBA as tuple of real numbers"
        return (self.red, self.green, self.blue, self.alphaDflt)

    @property
    def rgb8(self):
        "RGB as tuple of 8-bit ints"
        return (_realToInt8(self.red), _realToInt8(self.green), _realToInt8(self.blue))

    @property
    def rgba8(self):
        "RGBA as tuple of 8-bit ints"
        return (_realToInt8(self.red), _realToInt8(self.green), _realToInt8(self.blue), _realToInt8(self.alphaDflt))

    @property
    def packRgb8(self):
        "return packed 8-bit int RGB values (e.g.#008000)"
        return (_realToInt8(self.red) << 16) | (_realToInt8(self.green) << 8) | _realToInt8(self.blue)

    @property
    def hsv(self):
        "HSV as tuple of real numbers"
        return (self.hue, self.saturation, self.value)

    @property
    def hsva(self):
        "HSVA as tuple of real numbers"
        return (self.hue, self.saturation, self.value, self.alphaDflt)

    @property
    def hsv8(self):
        "HSV as tuple of 8-bit integers"
        return (self.hue8, self.saturation8, self.value8)

    @property
    def hsva8(self):
        "HSVA as tuple of 8-bit integers"
        return (self.hue8, self.saturation8, self.value8, self.alpha8)

    def toHtmlColor(self):
        "as an html color"
        return "#%02x%02x%02x" % self.rgb8

    def toRgbStr(self, sep=",", pos=4):
        "convert to a string of real RGB values, separated by sep"
        return "%0.*f%s%0.*f%s%0.*f" % (pos, self.red, sep, pos, self.green, sep, pos, self.blue)

    def toRgbaStr(self, sep=",", pos=4):
        "convert to a string of real RGBA values, separated by sep"
        return "%s%s%0.*f" % (self.toRgbStr(sep, pos), sep, pos, self.alphaDflt)

    def toRgb8Str(self, sep=","):
        "convert to a string of 8-bit RGB values, separated by sep"
        return "{}{}{}{}{}".format(_realToInt8(self.red), sep, _realToInt8(self.green), sep, _realToInt8(self.blue))

    def toRgba8Str(self, sep=","):
        "convert to a string of 8-bit RGBA values, separated by sep"
        return "{}{}{}".format(self.toRgb8Str(sep), sep, _realToInt8(self.alphaDflt))

    def toHsvStr(self, sep=",", pos=4):
        "convert to a string of real HSV values, separated by sep"
        return "%0.*f%s%0.*f%s%0.*f" % (pos, self.hue, sep, pos, self.saturation, sep, pos, self.value)

    def toHsv8Str(self, sep=","):
        "convert to a string of integer HSV values, separated by sep"
        return "{}{}{}{}{}".format(self.hue8, sep, self.saturation8, sep, self.value8)

    def setRed(self, red):
        "Create a new Color object with red set to the specified real number"
        return self.fromRgb(red, self.green, self.blue, self.alpha)

    def setGreen(self, green):
        "Create a new Color object with green set to the specified real number"
        return self.fromRgb(self.red, green, self.blue, self.alpha)

    def setBlue(self, blue):
        "Create a new Color object with blue set to the new real number"
        return self.fromRgb(self.red, self.green, blue, self.alpha)

    def setAlpha(self, alpha):
        "Create a new Color object with blue set to the new real number"
        return self.fromRgb(self.red, self.green, self.blue, alpha)

    def setHue(self, hue):
        "Create a new Color object with hue set to the specified real number"
        return self.fromHsv(hue, self.saturation, self.value, self.alpha)

    def setSaturation(self, sat):
        "Create a new Color object with saturation set to the specified real number"
        return self.fromHsv(self.hue, sat, self.value, self.alpha)

    def setValue(self, val):
        "Create a new Color object with value set to the new real number"
        return self.fromHsv(self.hue, self.saturation, val, self.alpha)

    @staticmethod
    def fromRgb(r, g, b, a=None):
        "construct from real RGB or RGBA values"
        assert (0.0 <= r <= 1.0)
        assert (0.0 <= g <= 1.0)
        assert (0.0 <= b <= 1.0)
        assert (a is None) or (0.0 <= a <= 1.0)
        return Color(r, g, b, a)

    @staticmethod
    def fromRgb8(r, g, b, a=None):
        "construct from 8-bit int RGB values"
        a8 = _int8ToReal(a) if a is not None else None
        return Color.fromRgb(_int8ToReal(r), _int8ToReal(g), _int8ToReal(b), a8)

    @staticmethod
    def fromPackRgb8(c):
        "construct from packed 8-bit int RGB values (e.g. #008000)"
        r = (c >> 16) & 0xff
        g = (c >> 8) & 0xff
        b = c & 0xff
        return Color.fromRgb(_int8ToReal(r), _int8ToReal(g), _int8ToReal(b))

    @staticmethod
    def fromRgb8Str(rgb8str):
        "construct from 8-bit int RGB values in a comma separate string (100,200,50)"
        try:
            parts = rgb8str.split(',')
            if len(parts) != 3:
                raise ValueError("expected three comma separated integers: {}".format(rgb8str))
            for i in range(3):
                parts[i] = int(parts[i])
                if not (0 <= parts[i] <= 255):
                    raise ValueError("color must be in the range 0..255: {}", parts[i])
            return Color.fromRgb8(*parts)
        except ValueError as ex:
            raise ValueError("invalid RGB8 color string: {}".format(rgb8str)) from ex

    @staticmethod
    def fromHsv(h, s, v, a=None):
        "construct from real HSV values"
        assert (0.0 <= h <= 1.0)
        assert (0.0 <= s <= 1.0)
        assert (0.0 <= v <= 1.0)
        assert (a is None) or (0.0 <= a <= 1.0)
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        return Color(r, g, b, a, hsv=(h, s, v))

    @staticmethod
    def fromHsv8(h, s, v):
        "construct from integer HSV values"
        return Color.fromHsv(h / 306.0, s / 100.0, v / 100.0)

    @staticmethod
    def fromHtmlColor(hcolor):
        "construct from a HTML color string in the form #804c66"
        mat = re.match("^#([0-9a-fA-F]{6})$", hcolor)
        if not mat:
            raise ValueError("invalid HTML color: {}".format(hcolor))
        return Color.fromPackRgb8(int(mat.group(1), base=16))

    def distance(self, other):
        "calculate euclidean distance between two colors"
        return math.sqrt(pow(self.red - other.red, 2) +
                         pow(self.green - other.green, 2) +
                         pow(self.blue - other.blue, 2))

    def complementary(self):
        "return complementary color to this one"
        return Color.fromRgb(1.0 - self.red, 1.0 - self.green, 1.0 - self.blue, self.alpha)

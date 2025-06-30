# Copyright 2006-2025 Mark Diekhans
# Copyright of original unknown (https://goodcode.io/articles/python-dict-object/)

# dictionary with keys as object files, based on:
# https://goodcode.io/articles/python-dict-object/
from collections import defaultdict

def _attributeError(name):
    raise AttributeError(f"No such attribute: `{name}'")

class ObjDict(dict):
    """Dict object where keys are field names.
    This is useful for JSON by doing:
       json.load(fh, object_pairs_hook=ObjDict)

    When inserting a dict, it must be explicitly converted to an ObjDict if
    desired.
    """
    __slots__ = ()

    def __getattr__(self, name):
        if name not in self:
            _attributeError(name)
        return self[name]

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name not in self:
            _attributeError(name)
        del self[name]


class DefaultObjDict(defaultdict):
    """defaultdict-based object where keys are field names.
    This is useful for JSON by doing:
       json.load(fh, object_pairs_hook=defaultObjDictJsonHook(default_factory))

    When inserting a dict, it must be explicitly converted to an DefaultObjDict or ObjDict if
    desired.
    """
    __slots__ = ()

    def __getattr__(self, name):
        return self[name]

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name not in self:
            _attributeError(name)
        del self[name]


def defaultObjDictJsonHook(default_factory=None):
    """"factory to create DefaultObjDict objects from JSON"""
    from functools import partial   # rarely used, so lazy load
    return partial(DefaultObjDict, default_factory)

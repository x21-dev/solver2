"""The module's ctypes declarations held to the C header.

solver2/_c.py is written by hand from bindings/c/include/solver2.h. A
struct with a field too few, or a parameter of another width, would
corrupt memory without a word, so this reads the header and compares:
every constant's value, every struct's fields and their types in order,
every callback's and every function's result and parameters; and that
neither side has a name the other lacks.

The header is in a plain style (one declaration per statement, no macros
in declarations), which is all this parser reads.
"""
import ctypes as C
import os
import re
import unittest

from _paths import ROOT
from solver2 import _c

HEADER = os.path.join(ROOT, "bindings", "c", "include", "solver2.h")

SCALARS = {"int": C.c_int, "double": C.c_double, "size_t": C.c_size_t, "uint64_t": C.c_uint64}


def ctype(text):
    """A C type of the header as the ctypes type _c.py should have for it."""
    words = text.replace("*", " * ").split()
    stars = words.count("*")
    base = " ".join(w for w in words if w not in ("*", "const"))
    if base == "char" and stars >= 1:
        t, stars = C.c_char_p, stars - 1
    elif base == "void" and stars >= 1:
        t, stars = C.c_void_p, stars - 1
    elif base in _c.HANDLES and stars >= 1:         # a handle is a number to the caller
        t, stars = C.c_void_p, stars - 1
    elif base == "unsigned char":
        t = C.c_ubyte
    elif base in SCALARS:
        t = SCALARS[base]
    elif base in _c.STRUCTS:
        t = _c.STRUCTS[base]
    elif base in _c.CALLBACKS:
        t = _c.CALLBACKS[base]
    else:
        raise AssertionError("solver2.h has a type this test doesn't know: %r" % text)
    for _ in range(stars):
        t = C.POINTER(t)
    return t


def parameters(text):
    """The types and names of a parameter list."""
    text = " ".join(text.split())
    if text in ("", "void"):
        return []
    out = []
    for part in text.split(","):
        kind, name = re.fullmatch(r"(.*?)(\w+)", part.strip()).groups()
        out.append((name, ctype(kind)))
    return out


class Header:
    def __init__(self):
        with open(HEADER, encoding="utf-8") as f:
            text = re.sub(r"/\*.*?\*/", " ", f.read(), flags=re.S)
        self.constants = {name: int(value) for name, value in re.findall(r"#define (SOLVER2_\w+) (\d+)", text)}
        for body in re.findall(r"\benum\s*\{(.*?)\}", text, flags=re.S):
            self.constants.update((name, int(value)) for name, value in re.findall(r"(SOLVER2_\w+)\s*=\s*(\d+)",
                                                                                   body))
        self.handles = re.findall(r"typedef struct (\w+) \1;", text)
        self.structs = {}
        for name, body, alias in re.findall(r"typedef struct (\w+)\s*\{(.*?)\}\s*(\w+);", text, flags=re.S):
            assert name == alias
            fields = []
            for statement in body.split(";"):
                if statement.strip():
                    kind, field = re.fullmatch(r"(.*?)(\w+)", " ".join(statement.split())).groups()
                    fields.append((field, ctype(kind)))
            self.structs[name] = fields
        self.callbacks = {name: (ctype(result), parameters(args)) for result, name, args in
                          re.findall(r"typedef (\w+) \(\*(\w+)\)\((.*?)\);", text, flags=re.S)}
        self.functions = {name: (ctype(result), parameters(args)) for result, name, args in
                          re.findall(r"^((?:const )?\w+\*?) (solver2_\w+)\((.*?)\);", text, flags=re.S | re.M)}


class Conformance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.header = Header()

    def test_the_header_was_read(self):
        h = self.header
        self.assertEqual((len(h.functions), len(h.structs), len(h.callbacks), len(h.handles), len(h.constants)),
                         (38, 14, 3, 4, 10))

    def test_constants(self):
        self.assertEqual(self.header.constants, _c.CONSTANTS)

    def test_handles(self):
        self.assertEqual(tuple(self.header.handles), _c.HANDLES)

    def test_structs(self):
        self.assertEqual(sorted(self.header.structs), sorted(_c.STRUCTS))
        for name, fields in self.header.structs.items():
            self.assertEqual(fields, list(_c.STRUCTS[name]._fields_), name)

    def test_callbacks(self):
        self.assertEqual(sorted(self.header.callbacks), sorted(_c.CALLBACKS))
        for name, (result, args) in self.header.callbacks.items():
            mine = _c.CALLBACKS[name]
            self.assertEqual((result, tuple(t for _, t in args)), (mine._restype_, mine._argtypes_), name)

    def test_functions(self):
        self.assertEqual(sorted(self.header.functions), sorted(_c.PROTOTYPES))
        for name, (result, args) in self.header.functions.items():
            self.assertEqual((result, [t for _, t in args]), (_c.PROTOTYPES[name][0], _c.PROTOTYPES[name][1]),
                             name)

    def test_the_library_is_the_modules(self):
        import solver2
        self.assertEqual(solver2._lib.solver2_abi(), self.header.constants["SOLVER2_ABI"])
        for name, (result, args) in _c.PROTOTYPES.items():
            function = getattr(solver2._lib, name)
            self.assertEqual((function.restype, function.argtypes), (result, args), name)


if __name__ == "__main__":
    unittest.main()

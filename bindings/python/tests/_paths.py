"""Where things are, for the tests: the repository's root, and the module on the path."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON = os.path.normpath(os.path.join(HERE, os.pardir))
ROOT = os.path.normpath(os.path.join(PYTHON, os.pardir, os.pardir))

for path in (os.path.join(PYTHON, "examples"), PYTHON):
    if path not in sys.path:
        sys.path.insert(0, path)

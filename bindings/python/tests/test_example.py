"""The same bits through Python: the Python worked example
(examples/beam.py) against the D program it mirrors (docs/example-class),
which runs the same problem natively, as a D class.

The two must print the same lines, the "bits" lines included: every
threshold, every elite's value and genes and the QD score as bit patterns.
The D program is built here with dub, so the comparison is always of the
code as it stands, and on one platform (runs differ between platforms).
"""
import os
import shutil
import subprocess
import unittest

import _paths
import beam

TWIN = os.path.join(_paths.ROOT, "docs", "example-class")


class TheWorkedExample(unittest.TestCase):
    def test_python_prints_what_d_prints(self):
        dub = shutil.which("dub")
        self.assertIsNotNone(dub, "dub is needed to build docs/example-class, which this compares with")
        built = subprocess.run([dub, "build", "--compiler=ldc2", "--build=release"], cwd=TWIN, capture_output=True,
                               text=True)
        self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
        program = os.path.join(TWIN, "solver2-example-class" + (".exe" if os.name == "nt" else ""))
        ran = subprocess.run([program], capture_output=True, text=True)
        self.assertEqual(ran.returncode, 0, ran.stdout + ran.stderr)
        native = ran.stdout.splitlines()

        mine = []
        beam.main(out=mine.append)

        self.assertEqual(mine, native)
        self.assertEqual(len(mine), 16)
        self.assertEqual(sum(line.startswith("bits:") for line in mine), 2)
        self.assertIn("the resumed run's archive is the uninterrupted run's, bit for bit", mine)


if __name__ == "__main__":
    unittest.main()

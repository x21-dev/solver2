"""Checks an unpacked release archive where it lies: the module loads the
library that came with it, and the worked example runs to its end.

    python check.py <the unpacked folder> [--expected <file>]

With --expected, the example must print that file's lines exactly: what
docs/example-class printed where the archive was built. That holds on the
same platform and C library; elsewhere runs may differ in their last bits,
so there the check is only that the example runs and its resumed run
equals its uninterrupted one.

It needs nothing of the repository, so that it can run where there is no
checkout and no D compiler.
"""
import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser(description="Checks an unpacked solver2 release.")
    parser.add_argument("folder")
    parser.add_argument("--expected")
    args = parser.parse_args()
    folder = os.path.abspath(args.folder)

    if os.environ.get("SOLVER2_LIBRARY"):
        raise SystemExit("SOLVER2_LIBRARY is set: the check must load the archive's own library")
    sys.path.insert(0, os.path.join(folder, "examples"))
    sys.path.insert(0, folder)
    import solver2
    import beam

    loaded = os.path.dirname(os.path.realpath(solver2.library_path))
    if loaded != os.path.realpath(os.path.join(folder, "solver2")):
        raise SystemExit("the module loaded %s, which is not the archive's" % solver2.library_path)
    if not solver2.checked:
        raise SystemExit("the archive's library was built without bounds checks")
    print("loaded %s: interface version %d, bounds checks on" % (solver2.library_path, solver2.abi))

    lines = []
    beam.main(out=lines.append)
    print("\n".join(lines))
    if "the resumed run's archive is the uninterrupted run's, bit for bit" not in lines or len(lines) != 16:
        raise SystemExit("the example didn't run to its end")

    if args.expected:
        with open(args.expected, encoding="utf-8") as f:
            expected = f.read().splitlines()
        if lines != expected:
            for mine, theirs in zip(lines, expected):
                if mine != theirs:
                    print("differs:\n  here:     %s\n  expected: %s" % (mine, theirs))
            raise SystemExit("the example printed other lines than where the archive was built")
        print("the same %d lines as where the archive was built" % len(lines))
    print("the archive works")


if __name__ == "__main__":
    main()

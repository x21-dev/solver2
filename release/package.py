"""Makes a release archive of solver2 for the platform it runs on: the
shared library, what it needs beside it, the Python module, the C header
and the worked example, in one folder that works where it is unpacked.

    python release/package.py --version 0.1.0 [--out dist]

The library must be built first, in bindings/c, with bounds checks (the
usual build): `dub build --compiler=ldc2 --build=release`.

The archive, solver2-<version>-<platform>-x86_64.zip (Windows) or .tar.gz
(Linux), unpacks to a folder of that name:

    solver2/        the Python module, with the library inside it: put the
                    folder that holds it on Python's path
    include/        solver2.h, for callers in C
    examples/       beam.py, the worked example
    README.md       how to use it, and what it needs
    LICENSE.md

On Linux the library uses LDC's two shared runtime libraries. They are
copied from the LDC that built it, under the names the library asks for,
next to it: its run path starts at its own folder, so it finds them there.
"""
import argparse
import ctypes
import os
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WINDOWS = os.name == "nt"
LIBRARY = "solver2.dll" if WINDOWS else "libsolver2.so"
#: What the library needs beside it on Linux, by the start of its name.
RUNTIME = ("libdruntime-ldc-shared.so", "libphobos2-ldc-shared.so")


def fail(message):
    raise SystemExit("package.py: " + message)


def runtime_libraries(library):
    """The LDC runtime libraries `library` is linked against, as (the name
    it asks for, the file the loader finds)."""
    listed = subprocess.run(["ldd", library], capture_output=True, text=True, check=True).stdout
    found = []
    for line in listed.splitlines():
        match = re.match(r"\s*(\S+) => (\S+)", line)
        if match and match.group(1).startswith(RUNTIME):
            if not os.path.isfile(match.group(2)):
                fail("%s needs %s, which the loader doesn't find" % (library, match.group(1)))
            found.append((match.group(1), match.group(2)))
    if len(found) != len(RUNTIME):
        fail("%s should use LDC's shared runtime (%s), and ldd shows %s. A library with the runtime linked in "
             "statically crashes on Linux: build it as bindings/c/dub.json says."
             % (library, ", ".join(RUNTIME), [name for name, _ in found] or "neither"))
    return found


def describe(library):
    """The library's interface version and whether it checks bounds. Both
    answer before solver2_init."""
    lib = ctypes.CDLL(library)
    return lib.solver2_abi(), bool(lib.solver2_checked())


def main():
    parser = argparse.ArgumentParser(description="Makes solver2's release archive for this platform.")
    parser.add_argument("--version", required=True, help="the release's version, as in its tag without the v")
    parser.add_argument("--out", default=os.path.join(ROOT, "dist"), help="where to put it (default: dist)")
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9A-Za-z][0-9A-Za-z.+-]*", args.version):
        fail("%r is no version" % args.version)

    built = os.path.join(ROOT, "bindings", "c", LIBRARY)
    if not os.path.isfile(built):
        fail("%s is not there. Build it in bindings/c: dub build --compiler=ldc2 --build=release" % built)
    abi, checked = describe(built)
    if not checked:
        fail("%s was built without bounds checks. A release has them: build the usual configuration." % built)

    name = "solver2-%s-%s-x86_64" % (args.version, "windows" if WINDOWS else "linux")
    stage = os.path.join(args.out, name)
    if os.path.exists(stage):
        shutil.rmtree(stage)
    module = os.path.join(stage, "solver2")
    os.makedirs(module)
    os.makedirs(os.path.join(stage, "include"))
    os.makedirs(os.path.join(stage, "examples"))

    source = os.path.join(ROOT, "bindings", "python", "solver2")
    for file in sorted(os.listdir(source)):
        if file.endswith(".py"):
            shutil.copy(os.path.join(source, file), module)
    shutil.copy(built, module)
    beside = [LIBRARY]
    if WINDOWS:
        # The import library, for a C program that links against the DLL.
        shutil.copy(os.path.join(ROOT, "bindings", "c", "solver2.lib"), module)
        beside.append("solver2.lib")
    else:
        for asked, file in runtime_libraries(built):
            shutil.copy(file, os.path.join(module, asked))      # follows a link: the file, under the asked name
            beside.append(asked)
    shutil.copy(os.path.join(ROOT, "bindings", "c", "include", "solver2.h"), os.path.join(stage, "include"))
    shutil.copy(os.path.join(ROOT, "bindings", "python", "examples", "beam.py"), os.path.join(stage, "examples"))
    shutil.copy(os.path.join(ROOT, "LICENSE.md"), stage)

    with open(os.path.join(ROOT, "release", "README.md"), encoding="utf-8") as f:
        readme = f.read()
    for key, value in (("{version}", args.version), ("{abi}", str(abi)), ("{library}", LIBRARY),
                       ("{beside}", ", ".join("`%s`" % b for b in beside))):
        readme = readme.replace(key, value)
    left = re.findall(r"\{[a-z]+\}", readme)
    if left:
        fail("release/README.md has placeholders this doesn't fill: %s" % left)
    with open(os.path.join(stage, "README.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write(readme)

    if WINDOWS:
        archive = stage + ".zip"
        with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
            for folder, _, files in os.walk(stage):
                for file in sorted(files):
                    path = os.path.join(folder, file)
                    z.write(path, os.path.relpath(path, args.out))
    else:
        archive = stage + ".tar.gz"
        with tarfile.open(archive, "w:gz") as t:
            t.add(stage, arcname=name)
    print(archive)


if __name__ == "__main__":
    sys.exit(main())

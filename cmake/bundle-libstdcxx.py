"""Make an auditwheel-repaired Linux wheel run on older glibc (#122).

conda-forge's RDKit and Qt need a newer libstdc++ than older distros ship, which
pushes the wheel up to manylinux_2_39. auditwheel won't bundle libstdc++ (it is
on the manylinux allow-list), so this does it: copy conda's libstdc++ and
libgcc_s into the wheel under private names, point every ELF at them, and retag
the wheel for the newest GLIBC symbol still needed.

usage: python bundle-libstdcxx.py WHEEL CONDA_LIB OUT_DIR
"""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

wheel, conda_lib, out_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
RENAMES = {"libstdc++.so.6": "libstdc++-penzene.so.6", "libgcc_s.so.1": "libgcc_s-penzene.so.1"}


def run(*args):
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"{' '.join(map(str, args))} failed:\n{r.stdout}{r.stderr}")
    return r.stdout


def is_elf(p):
    with open(p, "rb") as f:
        return f.read(4) == b"\x7fELF"


with tempfile.TemporaryDirectory() as tmp:
    run(sys.executable, "-m", "wheel", "unpack", "-d", tmp, str(wheel))
    root = next(Path(tmp).iterdir())
    libs = next(root.glob("*.libs"))  # auditwheel's penzene.libs

    for old, new in RENAMES.items():
        shutil.copy(conda_lib / old, libs / new)  # follows the symlink to the real file
        run("patchelf", "--set-soname", new, str(libs / new))

    elves = [p for p in root.rglob("*") if p.is_file() and is_elf(p)]
    for p in elves:
        needed = run("patchelf", "--print-needed", str(p)).split()
        for old, new in RENAMES.items():
            if old in needed:
                run("patchelf", "--replace-needed", old, new, str(p))
        if p.parent == libs:  # bundled libs find each other next to themselves
            run("patchelf", "--set-rpath", "$ORIGIN", str(p))

    # Newest GLIBC symbol version anything still needs from the system.
    newest = (2, 17)
    for p in elves:
        for major, minor in re.findall(r"GLIBC_(\d+)\.(\d+)", run("objdump", "-T", str(p))):
            newest = max(newest, (int(major), int(minor)))
    tag = f"manylinux_{newest[0]}_{newest[1]}_x86_64"

    packed = Path(tmp) / "packed"
    packed.mkdir()
    run(sys.executable, "-m", "wheel", "pack", "-d", str(packed), str(root))
    built = next(packed.glob("*.whl"))
    run(sys.executable, "-m", "wheel", "tags", "--remove", "--platform-tag", tag, str(built))
    out_dir.mkdir(parents=True, exist_ok=True)
    for w in packed.glob("*.whl"):
        shutil.move(str(w), out_dir / w.name)
        print(out_dir / w.name)

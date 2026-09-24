"""Fuzz the drawing engine through the Python module.

    python tests/fuzz.py [--seed N] [--trials N] [--samples DIR]

Random ChemDraw hotkey sequences, SMILES round trips (.penz and MOL),
abbreviations, and optionally a folder of ChemDraw files. Anything that raises,
produces invalid chemistry from valid input, or lays bonds out badly after Clean
is reported; the exit status is non-zero if a crash-class problem was found.
Uses an installed `penzene`, or the build tree's (build/python) if present.
"""
import argparse
import glob
import os
import random
import sys
import tempfile
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build", "python"))
import penzene as pz  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--trials", type=int, default=3000)
ap.add_argument("--samples", help="folder of .cdxml files to open and render")
args = ap.parse_args()
random.seed(args.seed)
problems = []  # (kind, detail, fatal)


def note(kind, detail, fatal):
    problems.append((kind, detail, fatal))


def bond_lengths(d):
    for b in d.bonds:
        p, q = d.atoms[b.a], d.atoms[b.b]
        yield ((p.x - q.x) ** 2 + (p.y - q.y) ** 2) ** 0.5


ATOM_KEYS = list("1234567890avuzkKcnowqspflCbihdBSLmePAONFEZMQHY+-")
BOND_KEYS = list("123whHazv45678dbyDBlcr90")
for _ in range(args.trials):
    d = pz.Document()
    hot = (d.add_atom("C"), -1)
    seq = ""
    for _ in range(random.randint(1, 12)):
        k = random.choice(BOND_KEYS if hot[1] >= 0 else ATOM_KEYS)
        seq += k
        try:
            hot = d.hotkeys(hot[0], k, bond=hot[1] if hot[1] >= 0 else None)
        except ValueError:
            break  # not a hotkey there: fine
        except Exception as e:  # noqa: BLE001
            note("hotkey raised", f"{seq}: {e!r}", True)
            break
    smi = d.to_smiles()
    try:
        if smi:
            pz.from_smiles(smi)
    except ValueError:
        note("SMILES does not re-parse", f"{seq} -> {smi}", True)
    try:
        d.to_svg()
        d.clean()
        d.to_svg()
    except Exception as e:  # noqa: BLE001
        note("render/clean raised", f"{seq}: {e!r}", True)
    if smi and any(not 10 < L < 19 for L in bond_lengths(d)):  # valid chemistry should clean evenly
        note("uneven bonds after Clean", seq, False)

out = tempfile.mkdtemp()
for s in ["CC(=O)Oc1ccccc1C(=O)O", "C[C@H](N)C(=O)O", "c1ccc2ccccc2c1", "C1CC1", "[Na+].[Cl-]", "C=C=C",
          "C#CC#C", "O=S(=O)(O)O", "C1CCCCCCCCCCC1", "c1ccncc1", "[NH4+]", "CC(C)(C)[Si](C)(C)OC",
          "F/C=C/F", "F/C=C\\F", "C[N+](C)(C)C", "B1OC(C)(C)C(C)(C)O1", "O", "c1cc2ccc3cccc4ccc(c1)c2c34"]:
    want = pz.from_smiles(s).to_smiles()
    for ext in ("penz", "mol"):
        d = pz.from_smiles(s)
        d.save(os.path.join(out, "x." + ext))
        got = pz.from_smiles(pz.read(os.path.join(out, "x." + ext)).to_smiles()).to_smiles()
        if got != want:
            note(f"{ext} round trip changes chemistry", f"{s} -> {got}", True)
    d = pz.from_smiles(s)
    d.clean()
    if pz.from_smiles(d.to_smiles()).to_smiles() != want:
        note("Clean changes chemistry", s, True)

for ab in ["Me", "Et", "iPr", "tBu", "Ph", "OMe", "NO2", "CF3", "CN", "CO2Me", "CO2Et", "CO2H", "CHO", "Ac",
           "OAc", "N3", "Boc", "Cbz", "Fmoc", "Bn", "Bz", "MgBr", "SO2Me", "Ts", "Ms", "Tf", "OTf", "OTs", "TMS",
           "TBS", "OTBS", "PMB", "Bpin", "nBu", "Pr", "Cy"]:
    d = pz.Document()
    c = d.add_atom("C")
    d.add_bond(c, d.add_atom(ab, x=14.4))
    if not d.to_smiles() or not d.formula:
        note("abbreviation has no chemistry", ab, True)
    d.clean()
    if len(d.atoms) != 2:
        note("Clean expanded an abbreviation", ab, True)

if args.samples:
    for f in sorted(glob.glob(os.path.join(args.samples, "*.cdxml"))):
        try:
            pz.read(f).to_svg()
        except Exception as e:  # noqa: BLE001
            note("CDXML sample failed", f"{os.path.basename(f)}: {e!r}", True)

for (kind, fatal), n in Counter((k, f) for k, _, f in problems).most_common():
    example = next(d for k, d, f in problems if k == kind)
    print(f"{'FAIL' if fatal else 'note'} {n:5d}  {kind}  (e.g. {example})")
fatal = sum(f for _, _, f in problems)
print(f"fuzz: seed {args.seed}, {args.trials} hotkey sequences, {fatal} fatal problems")
sys.exit(1 if fatal else 0)

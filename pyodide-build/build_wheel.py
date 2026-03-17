"""Build a Pyodide wheel from the /out/rdkit directory."""
import os
import re
import zipfile

# Read RDKit version from CMakeLists.txt
with open("/src/rdkit/CMakeLists.txt") as f:
    txt = f.read()
m = re.search(
    r'set\(RDKit_Year "(\d+)"\).*set\(RDKit_Month "(\d+)"\).*set\(RDKit_Revision "(\d+)"\)',
    txt,
    re.S,
)
if not m:
    raise RuntimeError("Could not parse RDKit version from CMakeLists.txt")
version = f"{m.group(1)}.{m.group(2)}.{m.group(3)}"
print(f"RDKit version: {version}")

whl_tag = "cp313-cp313-pyodide_2025_0_wasm32"
whl_name = f"rdkit_pyodide-{version}-{whl_tag}"
dist_info = f"/out/{whl_name}.dist-info"

# Create dist-info
os.makedirs(dist_info, exist_ok=True)

with open(f"{dist_info}/METADATA", "w") as f:
    f.write(
        f"Metadata-Version: 2.1\n"
        f"Name: rdkit-pyodide\n"
        f"Version: {version}\n"
        f"Summary: RDKit cheminformatics library for Pyodide\n"
    )

with open(f"{dist_info}/WHEEL", "w") as f:
    f.write(
        f"Wheel-Version: 1.0\n"
        f"Generator: rdkit-pyodide-build\n"
        f"Root-Is-Purelib: false\n"
        f"Tag: {whl_tag}\n"
    )

with open(f"{dist_info}/top_level.txt", "w") as f:
    f.write("rdkit\n")

# Build RECORD (empty — wheel spec allows omitting hashes for non-signed wheels)
with open(f"{dist_info}/RECORD", "w") as f:
    pass

# Create wheel zip
whl_path = f"/out/{whl_name}.whl"
with zipfile.ZipFile(whl_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for top_dir in ["rdkit", f"{whl_name}.dist-info"]:
        base = f"/out/{top_dir}"
        for root, _, files in os.walk(base):
            for fname in files:
                full = os.path.join(root, fname)
                arcname = os.path.relpath(full, "/out")
                zf.write(full, arcname)

size_mb = os.path.getsize(whl_path) / 1024 / 1024
print(f"Wheel: {whl_path} ({size_mb:.1f} MB)")

"""Snap smoke test: run the snap's own Python backend inside the snap
runtime and exercise the JSON IPC protocol exactly as the Tauri bridge
does. Verifies engines, staged python deps, and backend responsiveness.
"""
import json
import os
import subprocess
import sys

snap = os.environ["SNAP"]
assert snap, "must run inside snap runtime"

# Path checks for the wired integration points.
checks = {
    "wrapper": snap + "/bin/ligora-wrapper",
    "gpu-wrapper": snap + "/bin/gpu-2604-wrapper",
    "desktop-launch fwd": snap + "/snap/command-chain/desktop-launch",
    "run fwd": snap + "/snap/command-chain/run",
    "configure hook": snap + "/snap/hooks/configure",
    "platform launcher": snap + "/gnome-platform/command-chain/desktop-launch",
    "backend": snap + "/usr/lib/ligora-backend/ligora_backend/server.py",
    "python": snap + "/usr/bin/python3",
    "plipcmd": snap + "/usr/bin/plipcmd",
    "vina": snap + "/usr/bin/vina",
    "obabel": snap + "/usr/bin/obabel",
    "gmx": snap + "/usr/bin/gmx",
    "dssp": snap + "/usr/bin/dssp",
    "gmxdata": snap + "/usr/share/gromacs",
}
missing = [k for k, v in checks.items() if not os.path.exists(v)]
print("PATH CHECKS:", "ALL OK" if not missing else f"MISSING: {missing}")

env = dict(os.environ)
env.setdefault("LIGORA_BASE_DIR", os.environ.get("SNAP_USER_COMMON", "/tmp"))
env["LIGORA_BACKEND_DIR"] = snap + "/usr/lib/ligora-backend"
env["LIGORA_PYTHON"] = snap + "/usr/bin/python3"
env["PLIP_PATH"] = snap + "/usr/bin/plipcmd"
env["GMXDATA"] = snap + "/usr/share/gromacs"
env["BABEL_DATADIR"] = snap + "/usr/share/openbabel"
# mirror ligora-wrapper: Open Babel format plugins dir
import glob as _glob
_babel = sorted(_glob.glob(snap + "/usr/lib/x86_64-linux-gnu/openbabel/*/"))
if _babel:
    env["BABEL_LIBDIR"] = _babel[-1]
env["LD_LIBRARY_PATH"] = ":".join(
    p for p in [
        snap + "/usr/lib",
        snap + "/usr/lib/x86_64-linux-gnu",
        snap + "/lib/x86_64-linux-gnu",
        os.environ.get("LD_LIBRARY_PATH", ""),
    ] if p
)
env["PYTHONPATH"] = env["LIGORA_BACKEND_DIR"]

err_log = open("/tmp/ligora-backend-smoke.err", "w")
proc = subprocess.Popen(
    [env["LIGORA_PYTHON"], "-u", "-m", "ligora_backend.server"],
    env=env, cwd=env["LIGORA_BACKEND_DIR"],
    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
    stderr=err_log, text=True,
)

def ask(cmd, timeout=90):
    proc.stdin.write(json.dumps(cmd) + "\n")
    proc.stdin.flush()
    line = proc.stdout.readline()
    if not line:
        return None
    return json.loads(line)

# 1) get_status: engine registry + live data-source health
resp = ask({"type": "get_status", "payload": {}, "session_id": "smoke",
            "command_id": "s1"}, timeout=120)
if resp is None or not resp.get("success"):
    print("get_status FAILED:", json.dumps(resp)[:400] if resp else "no reply")
    err_log.flush()
    print("backend stderr:", open("/tmp/ligora-backend-smoke.err").read()[-1500:])
    sys.exit(1)
d = resp.get("data", {})
engines = {k: v.get("available") for k, v in (d.get("engines") or {}).items()}
sources = d.get("data_sources") or {}
print("ENGINES:", engines)
print("DATA SOURCES:", json.dumps(sources)[:300])

# 2) a quick open_pdb_id round trip with a tiny real structure
resp2 = ask({"type": "open_pdb_id", "payload": {"pdb_id": "1CRN"},
             "session_id": "smoke", "command_id": "s2"}, timeout=180)
ok2 = bool(resp2 and resp2.get("success"))
name = ((resp2 or {}).get("data") or {}).get("structure", {}).get("name", "?")
print("open_pdb_id 1CRN:", "OK" if ok2 else f"FAILED {json.dumps(resp2)[:200]}", name)

if not ok2:
    proc.terminate()
    print("SMOKE: FAIL (open_pdb_id)")
    sys.exit(2)

# 3) open a metalloprotein with real ligands (1HRC: cytochrome c + heme)
resp3 = ask({"type": "open_pdb_id", "payload": {"pdb_id": "1HRC"},
             "session_id": "smoke", "command_id": "s3"}, timeout=180)
ok3 = bool(resp3 and resp3.get("success"))
ligands = [l for l in ((resp3 or {}).get("data", {}).get("structure", {})
                       .get("ligands") or [])
           if l.get("residue_name") not in ("HOH", "DOD", "WAT")]
print("open_pdb_id 1HRC:", "OK" if ok3 else "FAILED",
      f"({len(ligands)} non-water ligands)")

ok3b = None
if ok3 and ligands:
    lid = ligands[0]["id"]
    r = ask({"type": "select_ligand", "payload": {"ligand_id": lid},
             "session_id": "smoke", "command_id": "s4"}, timeout=120)
    ok_sel = bool(r and r.get("success"))
    print(f"select_ligand {lid}:",
          "OK" if ok_sel else f"FAILED {json.dumps(r)[:200]}")
    if ok_sel:
        r2 = ask({"type": "run_contact_analysis", "payload": {},
                  "session_id": "smoke", "command_id": "s5"}, timeout=300)
        ok3b = bool(r2 and r2.get("success"))
        contacts = ((r2 or {}).get("data") or {}).get("contacts", [])
        print("run_contact_analysis (PLIP in snap):",
              f"OK ({len(contacts)} contacts)" if ok3b
              else f"FAILED {json.dumps(r2)[:300]}")
elif ok3:
    ok3b = False
    print("no non-water ligands found in 1HRC?!")

proc.terminate()
required = [e for e in ("vina", "plip", "openbabel", "rdkit", "gromacs")
            if e in engines]
verdict = ((not missing) and required
           and all(engines.get(e) for e in required) and ok2 and ok3)
print("SMOKE:", "PASS" if verdict else "ATTENTION NEEDED")
sys.exit(0 if verdict else 2)

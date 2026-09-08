"""
Short MD / minimization runs via the real GROMACS engine.

GROMACS (gmx) is invoked as an external binary — never reimplemented, and
never simulated when absent. What this adapter supports, honestly labeled:

- `em`:  steepest-descent energy minimization (gmx mdrun with an .mdp
         integrator = steep) in vacuum, no solvation.
- `md`:  a short vacuum molecular-dynamics snippet (leap-frog MD with
         velocity generation, no thermostat coupling unless the user
         configures one).

This is an analysis aid (geometry relaxation, quick stability probe), NOT
production dynamics and NOT solvated/periodic simulation. Every parameter
(em tolerance, step count, timestep, force field) is a user configuration
value with GROMACS-documented defaults — none are chemistry rules invented
here. All outputs (topology, energies, trajectory, final coordinates) are
GROMACS's own files, kept in the session workspace.
"""

from __future__ import annotations

import os
import re
import subprocess
import textwrap
import time
from pathlib import Path
from typing import Any, Dict, Optional

from ..config import get_config
from ..schemas import Structure
from .pockets import _write_protein_pdb

# mdp templates. Values are filled from user/config parameters; the
# structure of each file is GROMACS's own mdp format.
_MDP_EM = textwrap.dedent("""\
    ; Energy minimization (steepest descent) — Ligora-generated mdp
    integrator               = steep
    emtol                    = {emtol}
    emstep                   = 0.01
    nsteps                   = {nsteps}
    nstlist                  = 10
    cutoff-scheme            = Verlet
    coulombtype              = cut-off
    rcoulomb                 = 1.0
    rvdw                     = 1.0
    pbc                      = xyz
    """)

_MDP_MD = textwrap.dedent("""\
    ; Short vacuum MD snippet — Ligora-generated mdp
    integrator               = md
    dt                       = {dt}
    nsteps                   = {nsteps}
    nstenergy                = {nstenergy}
    nstxout                  = {nstxout}
    nstlist                  = 10
    cutoff-scheme            = Verlet
    coulombtype              = cut-off
    rcoulomb                 = 1.0
    rvdw                     = 1.0
    pbc                      = xyz
    gen_vel                  = yes
    gen_temp                 = {gen_temp}
    gen_seed                 = {gen_seed}
    tcoupl                   = no
    pcoupl                   = no
    """)


# --- X11-probe guard ---------------------------------------------------
#
# mdrun's hardware detection loads libXNVCtrl/NVML, which enumerates local
# X displays directly (not via DISPLAY) by connect()ing to the abstract
# UNIX sockets /tmp/.X11-unix/X<n>. On desktops where a stale Xwayland
# listener exists (observed live: gnome-shell's @X1 with a full accept
# queue), that connect() blocks forever and mdrun never starts — even
# though this GROMACS build has GPU support disabled and would never use
# the display. The guard is a tiny LD_PRELOAD shim that makes connect()
# to a *wedged* X socket fail fast with ECONNREFUSED so the enumeration
# moves on; it changes nothing else. Only sockets that provably don't
# answer (connect would hang) are affected — verified live: with the shim,
# the same mdrun converges normally.
_X_GUARD_C = textwrap.dedent("""\
#define _GNU_SOURCE
#include <dlfcn.h>
#include <sys/socket.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <poll.h>
#include <fcntl.h>
#include <unistd.h>
    typedef int (*connect_fn)(int, const struct sockaddr *, socklen_t);
    static connect_fn real_connect = NULL;
    static int is_xsock(const struct sockaddr *addr, socklen_t len) {
        if (!addr || addr->sa_family != AF_UNIX) return 0;
        const char *p = addr->sa_data;
        if (p[0] != '\\0') return 0;
        const char *target = "/tmp/.X11-unix/X";
        const size_t tlen = strlen(target);
        if ((size_t)len < 3 + tlen) return 0;
        return strncmp(p + 1, target, tlen) == 0;
    }
    int connect(int fd, const struct sockaddr *addr, socklen_t len) {
        if (!real_connect) real_connect = dlsym(RTLD_NEXT, "connect");
        if (is_xsock(addr, len)) {
            /* Non-blocking probe: refuse fast when the listener is wedged
             * (accept queue full), otherwise do the real connect. */
            int flags = fcntl(fd, F_GETFL, 0);
            int original = -1;
            if (flags != -1) {
                original = flags & O_NONBLOCK;
                fcntl(fd, F_SETFL, flags | O_NONBLOCK);
            }
            int rc = real_connect(fd, addr, len);
            int err = errno;
            if (original != -1) fcntl(fd, F_SETFL, original);
            if (rc == 0) return 0;
            if (err == EINPROGRESS) {
                struct pollfd pfd = {.fd = fd, .events = POLLOUT};
                int pr = poll(&pfd, 1, 3000);
                if (pr == 0) { errno = ECONNREFUSED; return -1; }
                if (pr > 0) {
                    int soerr = 0; socklen_t slen = sizeof(soerr);
                    getsockopt(fd, SOL_SOCKET, SO_ERROR, &soerr, &slen);
                    if (soerr == 0) return 0;
                    errno = soerr;
                    return -1;
                }
                errno = ECONNREFUSED;
                return -1;
            }
            return rc;
        }
        return real_connect(fd, addr, len);
    }
    int __connect(int fd, const struct sockaddr *addr, socklen_t len) {
        return connect(fd, addr, len);
    }
""")


def _x_guard_so(workdir: Path) -> Optional[Path]:
    """Compile the X-probe guard shim once; None when no compiler exists."""
    import tempfile
    cache_dir = Path(tempfile.gettempdir()) / "ligora-gmx-guard"
    so = cache_dir / "no_wedged_x.so"
    if so.exists():
        return so
    src = cache_dir / "no_wedged_x.c"
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        src.write_text(_X_GUARD_C, encoding="utf-8")
        import subprocess as _sp
        proc = _sp.run(
            ["cc", "-shared", "-fPIC", "-O2", "-o", str(so), str(src)],
            capture_output=True, timeout=60)
        if proc.returncode == 0 and so.exists():
            return so
    except (OSError, _sp.SubprocessError):
        pass
    try:
        src.unlink(missing_ok=True)
    except OSError:
        pass
    return None


class MDAdapter:
    """Run GROMACS minimization / short-MD snippets on the structure."""

    def __init__(self, binary_path: Optional[str] = None):
        self.config = get_config()
        # An explicitly configured binary is honored as-is: when it does
        # not exist we must report that honestly, never silently swap in
        # whatever else happens to be on PATH.
        self._binary_override = binary_path or os.environ.get(
            "LIGORA_GMX_PATH")

    # ------------------------------------------------------------------
    # Availability
    # ------------------------------------------------------------------

    def _gmx(self) -> Optional[str]:
        from shutil import which
        if self._binary_override:
            return self._binary_override \
                if Path(self._binary_override).exists() else None
        if self.config.gromacs_executable:
            configured = self.config.gromacs_executable
            return configured if Path(configured).exists() else None
        for name in ("gmx", "gmx_mpi"):
            found = which(name)
            if found:
                return found
        return None

    def is_available(self) -> bool:
        gmx = self._gmx()
        if not gmx:
            return False
        try:
            proc = subprocess.run(
                [gmx, "--version"], capture_output=True, timeout=30,
                text=True)
            return proc.returncode == 0 or bool(proc.stdout)
        except (OSError, subprocess.TimeoutExpired):
            return False

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self, structure: Structure, workdir: Path,
            mode: str = "em",
            force_field: Optional[str] = None,
            water_model: Optional[str] = None,
            nsteps: Optional[int] = None,
            emtol: Optional[float] = None,
            dt: Optional[float] = None,
            gen_temp: Optional[float] = None,
            gen_seed: Optional[int] = None,
            progress_cb=None,
            timeout_seconds: Optional[int] = None,
            ) -> Dict[str, Any]:
        """
        Prepare inputs (pdb2gmx), run (grompp + mdrun) and report energies.

        Every step's failure is surfaced as an error string; results are
        never simulated. Returns a dict with GROMACS's own outputs.
        """
        gmx = self._gmx()
        if not gmx:
            return {
                "available": False,
                "error": ("GROMACS (gmx) is not installed. Install it "
                          "(e.g. 'apt install gromacs') or set "
                          "LIGORA_GMX_PATH — MD/minimization is never "
                          "simulated."),
            }

        mode = (mode or "em").lower()
        if mode not in ("em", "md"):
            return {"available": False,
                    "error": f"Unknown MD mode: {mode!r} (use 'em' or 'md')"}

        cfg = self.config
        force_field = force_field or cfg.md_force_field
        water_model = water_model or cfg.md_water_model
        timeout = timeout_seconds or cfg.md_timeout

        workdir.mkdir(parents=True, exist_ok=True)
        started = time.time()
        log: list = []

        # Environment for GROMACS stages: headless (no display), and with
        # the X-probe guard when it could be built — mdrun's hardware
        # detection otherwise hangs forever on a wedged X display socket
        # on some desktops (see _X_GUARD_C above).
        stage_env = {k: v for k, v in os.environ.items()
                     if k != "DISPLAY"}
        guard = _x_guard_so(workdir)
        if guard:
            stage_env["LD_PRELOAD"] = str(guard)

        def run_stage(cmd: list, stdin_text: Optional[str] = None,
                      label: str = "") -> subprocess.CompletedProcess:
            if progress_cb:
                progress_cb(label)
            log.append(f"$ {' '.join(cmd)}")
            # Run inside the workspace subdirectory: GROMACS writes its
            # relative-named outputs (processed.gro, topol.top, ...) next
            # to the inputs, and they stay in the session workspace.
            proc = subprocess.run(
                cmd,
                input=stdin_text,
                stdin=None if stdin_text is not None
                else subprocess.DEVNULL,
                capture_output=True, text=True,
                cwd=str(workdir), timeout=timeout, env=stage_env)
            if proc.stdout:
                log.append(proc.stdout[-4000:])
            if proc.stderr:
                log.append(proc.stderr[-2000:])
            if proc.returncode != 0:
                tail = (proc.stderr or proc.stdout or "")[-500:]
                raise RuntimeError(
                    f"GROMACS stage '{label}' failed (exit "
                    f"{proc.returncode}): {tail}")
            return proc

        try:
            # 1. Polymer-only PDB for pdb2gmx (same strict writer PLIP uses).
            protein_pdb = workdir / "protein_input.pdb"
            n_atoms = _write_protein_pdb(structure, protein_pdb)
            if n_atoms == 0:
                return {
                    "available": True,
                    "error": ("The structure has no polymer atoms; there "
                              "is nothing to minimize or simulate."),
                }

            # 2. Topology generation (GROMACS's own pdb2gmx).
            # Force field: an explicit user choice is passed verbatim (and a
            # wrong name fails honestly). When unset, GROMACS picks from its
            # own shipped force-field list — its engine-owned default,
            # selected with stdin "1"; the actual choice is read back from
            # the generated topology and recorded.
            pdb2gmx_cmd = [
                gmx, "pdb2gmx", "-f", protein_pdb.name,
                "-o", "processed.gro", "-p", "topol.top",
                "-water", water_model, "-ignh"]
            pdb2gmx_stdin = None
            if force_field:
                pdb2gmx_cmd += ["-ff", force_field]
            else:
                pdb2gmx_stdin = "1\n"
            run_stage(pdb2gmx_cmd, stdin_text=pdb2gmx_stdin,
                      label="pdb2gmx (topology)")

            chosen_ff = _force_field_from_topology(
                workdir / "topol.top") or force_field

            # 2b. Put the molecule in a large cubic box (1 nm padding).
            # Modern GROMACS requires periodic boundaries with Verlet
            # lists; a big isolated box is the standard vacuum setup —
            # no solvent, no added ions, nothing simulated.
            run_stage(
                [gmx, "editconf", "-f", "processed.gro",
                 "-o", "boxed.gro", "-d", "1.0", "-bt", "cubic"],
                label="editconf (vacuum box)")

            # 3. mdp from user/config parameters.
            if mode == "em":
                mdp = _MDP_EM.format(
                    emtol=float(emtol if emtol is not None else 1000.0),
                    nsteps=int(nsteps if nsteps is not None else 50000))
            else:
                mdp = _MDP_MD.format(
                    dt=float(dt if dt is not None else 0.001),
                    nsteps=int(nsteps if nsteps is not None else 5000),
                    nstenergy=100, nstxout=100,
                    gen_temp=float(gen_temp if gen_temp is not None
                                   else 298.0),
                    gen_seed=int(gen_seed if gen_seed is not None
                                 else int(time.time()) % 32768))
            mdp_path = workdir / f"{mode}.mdp"
            mdp_path.write_text(mdp, encoding="utf-8")

            # 4. Preprocess into a run input (tpr).
            run_stage(
                [gmx, "grompp", "-f", mdp_path.name,
                 "-c", "boxed.gro",
                 "-p", "topol.top",
                 "-o", f"{mode}.tpr"],
                label="grompp (run input)")

            # 5. Run the engine.
            run_stage(
                [gmx, "mdrun", "-deffnm", mode,
                 "-s", f"{mode}.tpr"],
                label=f"mdrun ({mode})")

            # 6. Extract energies with gmx energy. For EM the meaningful
            # series is Potential; for MD the thermodynamic set (Potential,
            # Kinetic, Total Energy, Temperature) is what a short stability
            # probe reports. gmx energy's interactive menu matches by number
            # or exact term; the menu spellings are its own ('Kinetic-En.',
            # 'Total-Energy', 'Temperature' — observed from its live menu),
            # so the selection uses those verbatim. 0 ends selection.
            edr = workdir / f"{mode}.edr"
            energies: Dict[str, Any] = {}
            if edr.exists():
                selected = (["Potential"] if mode == "em" else
                            ["Potential", "Kinetic-En.", "Total-Energy",
                             "Temperature"])
                selection = "\n".join(selected) + "\n0\n"
                run_stage(
                    [gmx, "energy", "-f", f"{mode}.edr",
                     "-o", "energy.xvg"],
                    stdin_text=selection,
                    label="gmx energy")
                energies = _parse_xvg(workdir / "energy.xvg")

            runtime = time.time() - started
            return {
                "available": True,
                "mode": mode,
                "engine": "GROMACS (external binary, not reimplemented)",
                "binary": gmx,
                "force_field": chosen_ff,
                "force_field_explicit": bool(force_field),
                "water_model": water_model,
                "nsteps": int(nsteps if nsteps is not None
                              else (50000 if mode == "em" else 5000)),
                "runtime_seconds": round(runtime, 2),
                "atom_count": n_atoms,
                "energies": energies,
                "files": {
                    "topology": str(workdir / "topol.top"),
                    "boxed_structure": str(workdir / "boxed.gro"),
                    "run_input": str(workdir / f"{mode}.tpr"),
                    "final_structure": str(workdir / f"{mode}.gro"),
                    "energies_edr": str(edr) if edr.exists() else None,
                    "energy_xvg": str(workdir / "energy.xvg")
                    if (workdir / "energy.xvg").exists() else None,
                    "trajectory": str(workdir / f"{mode}.trr")
                    if (workdir / f"{mode}.trr").exists() else None,
                    "log": str(workdir / f"{mode}.log")
                    if (workdir / f"{mode}.log").exists() else None,
                },
                "scope_note": (
                    "Vacuum minimization/MD snippet (isolated box, no "
                    "solvent) for analysis context — not solvated "
                    "production dynamics."
                    if mode == "md" else
                    "Vacuum energy minimization (isolated box, no solvent) "
                    "for analysis context."),
                "log": "\n".join(log),
            }
        except (RuntimeError, subprocess.TimeoutExpired, OSError) as e:
            return {
                "available": True,
                "mode": mode,
                "error": str(e),
                "runtime_seconds": round(time.time() - started, 2),
                "log": "\n".join(log),
            }


def _force_field_from_topology(path: Path) -> Optional[str]:
    """
    Read the force field GROMACS actually used from its own generated
    topology (the forcefield.itp include names it). Real engine output.
    """
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    match = re.search(r'"(\S+?)\.ff[/\\]forcefield\.itp"', text)
    return match.group(1) if match else None


def _parse_xvg(path: Path) -> Dict[str, Any]:
    """
    Parse GROMACS's energy.xvg: the final row of real numeric series.
    Returns the last time and the selected series values by name.
    """
    names: list = []
    last_values: list = []
    last_time: Optional[float] = None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith("@ s") and "legend" in line:
            names.append(line.rsplit("legend ", 1)[-1].strip('"'))
        elif line.startswith(("#", "@")) or not line:
            continue
        else:
            parts = line.split()
            try:
                last_time = float(parts[0])
                last_values = [float(v) for v in parts[1:]]
            except ValueError:
                continue
    series = {n: v for n, v in zip(names, last_values) if names}
    return {"time": last_time, "series": series}

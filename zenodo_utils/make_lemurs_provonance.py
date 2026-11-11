#!/usr/bin/env python3
"""Build the LEMURS provenance package: random-engine states + merge scripts.

For each of the five LEMURS productions it collects
    <BASE>/<production>/randomSeeds/<CLUSTER>.<PROC>/currentRun.rndm
    <BASE>/<production>/merge.sh
reads merge.sh (with `hadd` replaced by a recorder) to find the jobs merged into
each released file, packages only their states, checks for duplicates, and writes
    LEMURS_provenance/README.md          what is inside and how to use it
    LEMURS_provenance/lemurs_seeds.csv   one row per job (state hash)
    LEMURS_provenance/MANIFEST.sha256    checksum of every packaged file
    LEMURS_provenance/<production>/...   the states and merge.sh, original layout
into one tarball. Standard library only; nothing on EOS is modified.

Usage:
    python3 make_lemurs_provenance.py                       # build the tarball
    python3 make_lemurs_provenance.py --dry-run             # checks + report only
    python3 make_lemurs_provenance.py --base /eos/... --out my.tar.gz
Exit codes: 0 ok, 1 duplicates or missing inputs (nothing written unless --force).
"""
import argparse, csv, hashlib, io, os, re, subprocess, sys, tarfile, tempfile, time
from collections import Counter, defaultdict
from pathlib import Path

BASE = "/eos/geant4/fastSim/ddfastsim"
PRODUCTIONS = {   # detector label -> production directory below BASE
    "Par04_SiW":     "/eos/user/p/praikwar/ddsim_out/Par04Flat",   # absolute: outside BASE
    "Par04_SciPb":   "Par04SciPb/dataset2_SciPb_1GeV1TeVFlat_theta0p87to2p27_phiFull",
    "ODD":           "ODD/dataset2_1GeV1TeVFlat_theta0p87to2p27_phiFull",
    "FCCee_CLD":     "FCCeeCLD/dataset2_1GeV100GeVFlat_theta0p87to2p27_phiFull",
    "FCCee_ALLEGRO": "FCCeeALLEGRO/dataset2_1GeV100GeVFlat_theta0p87to2p27_phiFull",
}
# Productions whose merge script was not kept: assume the procedure of the other Par04
# productions (10 parts x 200 jobs, ascending job order) and ship a reconstructed script.
ASSUMED_MERGE = {
    "Par04_SiW": {"cluster": "12388001",
                  "stem": "ddsim_mesh_Par04_gamma_500events_1GeV1TeV_GPSflat_edm4hep",
                  "parts": 10, "jobs_per_part": 200},
}
RECONSTRUCTED = """# merge_reconstructed.sh -- RECONSTRUCTED, NOT THE ORIGINAL SCRIPT
# The original merge script of this production was not kept. This file assumes the
# same procedure as the other Par04 productions (see Par04SciPb/.../merge.sh):
# {parts} parts of {jobs_per_part} jobs each, merged with hadd in ascending job order.
for i in {{1..{parts}}}; do
    start=$(({jobs_per_part} * (i - 1)))
    end=$(({jobs_per_part} * i - 1))
    eval "files={stem}_{cluster}.{{$start..$end}}.root"
    eval "hadd {stem}_{cluster}_part$i.root $files" &
done
wait
"""
DATASET_DOI = "10.5281/zenodo.17045562"
SOFTWARE_DOI = "10.5281/zenodo.17076360"
VERSIONS = ("key4hep release 2023-11-23 (build 2023-11-30-gyuooo, x86_64-almalinux9-gcc11.3.1-opt): "
            "Geant4 11.1.2, DD4hep 1.27.1, EDM4hep 0.10.2, ROOT 6.28.06; ddfastsim v1.0.0")
JOB_RE = re.compile(r"^(\d+)\.(\d+)$")
OUTPUT_RE = re.compile(r"_(\d+)\.(\d+)\.root$")        # per-job ddsim output, not merged parts
MIXMAX_RE = re.compile(r"V\[N\]\s*=\s*\{([^}]*)\}\s*;\s*counter\s*=\s*(\d+)", re.S)


def state_hash(text):
    m = MIXMAX_RE.search(text)
    if m:
        vals = tuple(int(v) for v in re.findall(r"-?\d+", m.group(1)))
        return "mixmax", hashlib.sha256(repr((vals, int(m.group(2)))).encode()).hexdigest()
    return "raw", hashlib.sha256(" ".join(text.split()).encode()).hexdigest()


def merged_jobs(pdir, merge):
    """Run merge.sh with `hadd` replaced by a recorder; return {output part: [input names]}.

    Only `hadd` is intercepted; the script's own echo/eval/brace expansion run as written,
    in the production directory, so the file lists are exactly what hadd received.
    """
    with tempfile.TemporaryDirectory() as tmp:
        stub = ('hadd() { out="$1"; shift; printf "%s\\n" "$@" > "$HADD_LOG/$(basename "$out").list"; }; '
                'export -f hadd 2>/dev/null; source ./merge.sh > /dev/null; wait')
        subprocess.run(["bash", "-c", stub], cwd=pdir, env={**os.environ, "HADD_LOG": tmp},
                       check=True, timeout=600)
        parts = {}
        for f in sorted(Path(tmp).glob("*.list")):
            parts[f.name[:-5]] = [l for l in f.read_text().split("\n") if l]
        return parts


def survey(label, pdir, allow_missing_merge=False, assumed=None):
    """States, merge lists and consistency of one production."""
    info = {"label": label, "dir": pdir, "states": [], "problems": [], "formats": Counter(),
            "parts": {}, "merged": [], "missing_inputs": [], "merged_without_state": [],
            "unmerged_states": 0}
    seeds, merge = pdir / "randomSeeds", pdir / "merge.sh"
    info["merge"] = merge if merge.is_file() else None
    if not seeds.is_dir():
        info["problems"].append(f"missing {seeds}")
    info["inferred"] = None
    info["assumed"] = None
    if info["merge"] is None and assumed:
        a = assumed
        info["assumed"] = a
        info["reconstructed_text"] = RECONSTRUCTED.format(**a)
        info["parts"] = {f"part{i}": [f"{a['stem']}_{a['cluster']}.{j}.root"
                                      for j in range(a["jobs_per_part"] * (i - 1), a["jobs_per_part"] * i)]
                         for i in range(1, a["parts"] + 1)}
    if info["merge"] is None and not allow_missing_merge and not assumed:
        info["problems"].append(f"missing {merge} (use --allow-missing-merge {label} to infer the job list)")
    if info["problems"]:
        return info
    if info["merge"] is None and not assumed:
        # no merge script kept: take every job of the largest cluster as merged
        jobs = [d.name for d in seeds.iterdir() if JOB_RE.match(d.name)]
        main = Counter(j.split(".")[0] for j in jobs).most_common(1)[0][0] if jobs else None
        sel = sorted((j for j in jobs if j.split(".")[0] == main), key=lambda j: int(j.split(".")[1]))
        info["parts"] = {"(inferred)": [f"inferred_{j}.root" for j in sel]}
        info["inferred"] = main
    elif info["merge"] is not None:
        try:
            info["parts"] = merged_jobs(pdir, merge)
        except Exception as e:                   # noqa: BLE001
            info["problems"].append(f"could not evaluate merge.sh: {e}")
            return info
    merged, info["jobs_per_part"] = [], {}
    for part, inputs in info["parts"].items():
        info["jobs_per_part"][part] = len(inputs) if (info["inferred"] or info["assumed"]) else \
            sum(1 for n in inputs if "*" not in n and (pdir / n).exists())
        for name in inputs:
            m = OUTPUT_RE.search(name)
            if not (info["inferred"] or info["assumed"]) and ("*" in name or not (pdir / name).exists()):
                info["missing_inputs"].append(f"{part}: {name}")
            if m and "*" not in name:
                merged.append(f"{m.group(1)}.{m.group(2)}")
    info["merged"] = merged
    states = {}
    for f in seeds.glob("*/currentRun.rndm"):
        if JOB_RE.match(f.parent.name):
            states[f.parent.name] = f
    info["unmerged_states"] = len(set(states) - set(merged))
    for job in merged:
        f = states.get(job)
        if f is None:
            info["merged_without_state"].append(job)
            continue
        text = f.read_text(errors="replace")
        if not text.strip():
            info["problems"].append(f"empty state file: {f}")
            continue
        fmt, h = state_hash(text)
        info["formats"][fmt] += 1
        info["states"].append((job, h, f))
    return info


def readme(infos, dup_groups, created):
    lines = [
        "# LEMURS provenance records",
        "",
        f"Random-engine states and merge scripts for the LEMURS dataset (doi:{DATASET_DOI}),",
        f"produced with ddfastsim (doi:{SOFTWARE_DOI}). Created {created} by make_lemurs_provenance.py.",
        "Licence: CC BY 4.0.",
        "",
        "## Contents",
        "",
        "- `<production>/randomSeeds/<CLUSTER>.<PROC>/currentRun.rndm`: state of the Geant4 MixMax",
        "  random-number engine at the start of each HTCondor job (saved with `/random/setSavingFlag 1`).",
        "  The seed was not fixed in the steering files; ddsim chose it for each job (standard DD4hep procedure).",
        "- `<production>/merge.sh`: the `hadd` commands that merged the per-job EDM4hep outputs into the",
        "  ten training files of each detector.",
        "- `lemurs_seeds.csv`: one row per job with the SHA-256 of its engine state.",
        "- `MANIFEST.sha256`: checksum of every file in this package.",
        "",
        "## Software",
        "",
        VERSIONS + ".",
        "",
        "## Reproducing a job",
        "",
        "Restore the saved state before the run starts (for example with the Geant4 command",
        "`/random/resetEngineFrom currentRun.rndm`) and run the job's ddfastsim configuration",
        "with the software above.",
        "",
        "## From a shower to its job",
        "",
        "Each `merge.sh` passes the job outputs to `hadd` in ascending order, 200 jobs of 500 events per",
        "training file: part *i* contains jobs 200(i-1) ... 200i-1. In a file without the event selection",
        "of `edm4hep2h5.py`, shower *k* of part *i* comes from job 200(i-1) + floor(k/500), event k mod 500.",
        "Files that passed the selection (FCCee_CLD parts 5-10, FCCee_ALLEGRO) keep this order, but the",
        "removed events must be accounted for by re-applying `edm4hep2h5.py` to the merged EDM4hep file.",
        "Check each `merge.sh` for productions with resubmitted jobs (see below).",
        "",
        "## Per production",
        "",
        "| Production | Jobs merged | Merged jobs with saved state | Merged jobs without state | Jobs per part | Inputs named in merge.sh but absent | Unmerged states (not packaged) |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for i in infos:
        per_part = ", ".join(f"{p.split('part')[-1].lstrip('_.').split('.')[0]}: {n}" for p, n in i.get("jobs_per_part", {}).items())
        lines.append(f"| {i['label']} | {len(i['merged'])} | {len(i['states'])} | "
                     + (", ".join(i["merged_without_state"]) or "none") + f" | {per_part} | "
                     + ("; ".join(i["missing_inputs"]) or "none") + f" | {i['unmerged_states']} |")
    for i in infos:
        if i.get("assumed"):
            a = i["assumed"]
            lines.append(f"\n{i['label']}: the original merge script was not kept. The package contains "
                         f"`merge_reconstructed.sh`, which assumes the procedure of the other Par04 productions "
                         f"({a['parts']} parts of {a['jobs_per_part']} jobs of cluster {a['cluster']}, ascending order); "
                         f"the shower-to-job mapping above holds under that assumption.")
        if i.get("inferred"):
            lines.append(f"\n{i['label']}: the merge script was not kept; the merged jobs are taken to be all "
                         f"jobs of cluster {i['inferred']} with a saved state.")
    lines += ["", "No two saved states of merged jobs are identical." if not dup_groups
              else f"WARNING: {len(dup_groups)} engine state(s) shared by several jobs; see lemurs_seeds.csv.", ""]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=BASE, type=Path, help=f"production base directory (default {BASE})")
    ap.add_argument("--out", default="LEMURS_provenance.tar.gz", type=Path, help="output tarball")
    ap.add_argument("--dry-run", action="store_true", help="run the checks and print the README only")
    ap.add_argument("--force", action="store_true", help="write the package even if checks fail")
    ap.add_argument("--allow-missing-merge", action="append", default=[], metavar="LABEL",
                    help="production without merge.sh: infer merged jobs from its largest cluster")
    ap.add_argument("--set", action="append", default=[], metavar="LABEL=RELPATH",
                    help="override a production directory, e.g. Par04_SiW=Par04/<dir>")
    args = ap.parse_args()

    if args.out.exists() and not args.dry_run:
        sys.exit(f"{args.out} exists; refusing to overwrite.")

    prods = dict(PRODUCTIONS)
    for kv in args.set:
        k, _, v = kv.partition('=')
        if k not in prods or not v:
            sys.exit(f'--set needs one of {list(prods)}=RELPATH')
        prods[k] = v
    infos = [survey(label, args.base / rel, label in args.allow_missing_merge, ASSUMED_MERGE.get(label))
             for label, rel in prods.items()]
    by_hash = defaultdict(list)
    for i in infos:
        for job, h, f in i["states"]:
            by_hash[h].append((i["label"], job))
    dup_groups = {h: v for h, v in by_hash.items() if len(v) > 1}

    for i in infos:
        print(f"{i['label']:14s} merged jobs {len(i['merged']):5d}  with state {len(i['states']):5d}  "
              f"without state {len(i['merged_without_state']):4d}  unmerged states {i['unmerged_states']:4d}  "
              f"formats {dict(i['formats'])}")
        if i["parts"]:
            print("   jobs per part: " + ", ".join(str(n) for n in i["jobs_per_part"].values()))
        for m in i["missing_inputs"]:
            print(f"   input named in merge.sh but absent: {m}")
        for p in i["problems"]:
            print(f"   PROBLEM: {p}")
    tot_m = sum(len(i["merged"]) for i in infos); tot_s = sum(len(i["states"]) for i in infos)
    print(f"TOTAL: {tot_s} of {tot_m} merged jobs have a saved state")
    print("Duplicates:", "none" if not dup_groups else f"{len(dup_groups)} shared state(s)")
    for h, v in dup_groups.items():
        print("   ", h[:12], v)

    failed = bool(dup_groups) or any(i["problems"] for i in infos)
    created = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())
    text = readme(infos, dup_groups, created)
    if args.dry_run:
        print("\n" + text)
        sys.exit(1 if failed else 0)
    if failed and not args.force:
        sys.exit("Checks failed; nothing written (use --force to package anyway).")

    # seeds CSV
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["production", "job", "state_sha256", "shared_by_jobs"])
    for i in infos:
        for job, h, _ in i["states"]:
            w.writerow([i["label"], job, h, len(by_hash[h])])
    generated = {"README.md": text.encode(), "lemurs_seeds.csv": buf.getvalue().encode()}

    # files to pack, keeping the original layout below BASE
    members = []
    for i in infos:
        try:
            rel = i["dir"].relative_to(args.base)
        except ValueError:                       # e.g. Par04_SiW in another EOS area
            rel = Path(i["label"], i["dir"].name)
        members += [(f, Path(rel, "randomSeeds", f.parent.name, f.name)) for _, _, f in i["states"]]
        if i["merge"]:
            members.append((i["merge"], Path(rel, "merge.sh")))
        elif i.get("reconstructed_text"):
            generated[str(Path(rel, "merge_reconstructed.sh"))] = i["reconstructed_text"].encode()
    manifest = [f"{hashlib.sha256(src.read_bytes()).hexdigest()}  {arc}" for src, arc in members]
    manifest += [f"{hashlib.sha256(data).hexdigest()}  {name}" for name, data in generated.items()]
    generated["MANIFEST.sha256"] = ("\n".join(manifest) + "\n").encode()

    root = "LEMURS_provenance"
    with tarfile.open(args.out, "w:gz") as tar:
        for name, data in generated.items():
            ti = tarfile.TarInfo(f"{root}/{name}"); ti.size = len(data); ti.mtime = int(time.time())
            tar.addfile(ti, io.BytesIO(data))
        for src, arc in members:
            tar.add(src, arcname=f"{root}/{arc}")
    print(f"\nWrote {args.out} ({args.out.stat().st_size / 1e6:.2f} MB, {len(members)} files + README, CSV, manifest)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Convert existing LEMURS FCC-ee barrel EDM4hep ROOT events to *_v2.h5.

Four original datasets only: showers (S,9,16,45), incident_energy, incident_theta,
incident_phi; float32, MeV/radians, gzip 9. Selection uses primary endpoint plus
EARLY saved-secondary vertices, NOT MCParticles.size()==1. See LEMURS_v2_README.md.
Streams TTree entries and bounded HDF5 batches; no RDataFrame/AsNumpy.
Hit copying is always scalar/indexed Python: no custom C++ compilation and no
ROOT-vector-to-NumPy conversion. PyROOT itself still uses runtime bindings.
A small .selection.npz stores scalar diagnostics for ALL original ROOT rows, so
plots can compare inclusive, numMC==1, and v2 on exactly the same original events.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import traceback

import h5py
import numpy as np

VERSION = 'lemurs-v2-converter-1.0.1'
GEOMETRY_COMMIT = 'a25ddc8d9bcfd5610e389cda487d3f00172e9443'
# Frozen from the uploaded k4geo-FastSim archive, NOT current upstream master.
# ALLEGRO barrel region includes the cryostat: 2100+49 = 2149 mm.
# ALLEGRO endcap region is EMEC: z=3200+50+10=3260, r=350+50+10=410.
GEOMETRIES = {
    'CLD': dict(name='FCCee_o2_v04_noTracker', sides=12,
                barrel_r=2150., barrel_outer_r=2352., barrel_half_z=2210.,
                endcap_z=2307., endcap_back_z=2509., endcap_inner_r=340.,
                endcap_outer_r=2350.),
    'ALLEGRO': dict(name='ALLEGRO_o1_v02_noTracker', sides=0,
                    barrel_r=2149., barrel_outer_r=2721., barrel_half_z=3100.,
                    endcap_z=3260., endcap_back_z=3710., endcap_inner_r=410.,
                    endcap_outer_r=2750.),
}
REASONS = {
    0: 'accepted',
    1: 'unexpected_generated_particles_or_origin',
    2: 'no_nominal_ECAL_intersection',
    3: 'primary_photon_endpoint_upstream',
    4: 'early_simulated_particle_produced_upstream',
}
C = 299.792458  # mm/ns
_NORMALS = [(math.cos(k*math.pi/6), math.sin(k*math.pi/6)) for k in range(12)]


def support(x, y, g):
    """CLD: maximum projection on the 12 face normals; ALLEGRO: radius."""
    if g['sides'] == 0:
        return math.hypot(x, y)
    return max(x*nx + y*ny for nx, ny in _NORMALS)


def upstream(x, y, z, g, margin):
    return (support(x, y, g) < g['barrel_r'] - margin
            and abs(z) < g['endcap_z'] - margin)


def entry_distance(px, py, pz, g):
    """First nominal straight IP ray crossing of barrel side/endcap FRONT.

    A central-hole crossing is NOT an endcap entry. Not a Geant4 navigator;
    this intentionally retains the already-tested CLD nominal-surface rule.
    """
    norm = math.sqrt(px*px + py*py + pz*pz)
    if not math.isfinite(norm) or norm <= 0:
        return math.inf
    dx, dy, dz = px/norm, py/norm, pz/norm
    first = math.inf
    transverse = support(dx, dy, g)
    if transverse > 0:
        sb = g['barrel_r']/transverse
        if abs(sb*dz) <= g['barrel_half_z']:
            first = sb
    if abs(dz) > 1e-12:
        se = g['endcap_z']/abs(dz)
        sr = support(se*dx, se*dy, g)
        if g['endcap_inner_r'] <= sr <= g['endcap_outer_r']:
            first = min(first, se)
    return first


def classify(particles, g, boundary_margin=1., time_margin=.1):
    """EDM4hep data-record objects; returns one mutually exclusive reason."""
    if len(particles) == 0:
        return 1
    p = particles[0]
    if int(p.PDG) != 22 or int(p.generatorStatus) != 1:
        return 1
    for v in (p.vertex, p.endpoint, p.momentum):
        if not all(math.isfinite(float(getattr(v, k))) for k in ('x','y','z')):
            raise ValueError('Nonfinite primary truth coordinates/momentum')
    if not math.isfinite(float(p.time)):
        raise ValueError('Nonfinite primary time')
    if max(abs(float(getattr(p.vertex, k))) for k in ('x','y','z')) > 1.:
        return 1
    for i in range(1, len(particles)):
        if int(particles[i].generatorStatus) != 0:
            return 1
    length = entry_distance(float(p.momentum.x), float(p.momentum.y),
                            float(p.momentum.z), g)
    if not math.isfinite(length):
        return 2
    if upstream(float(p.endpoint.x), float(p.endpoint.y), float(p.endpoint.z),
                g, boundary_margin):
        return 3
    entry_time = float(p.time) + length/C
    for i in range(1, len(particles)):
        s = particles[i]
        vals = [float(s.vertex.x), float(s.vertex.y), float(s.vertex.z), float(s.time)]
        if not all(math.isfinite(x) for x in vals):
            raise ValueError('Nonfinite saved-secondary vertex/time')
        if vals[3] < entry_time - time_margin and upstream(*vals[:3], g, boundary_margin):
            return 4
    return 0


def project_hits(hits, shape, sampling_fraction):
    """ROOT hit positions encode (r,phi,z) BIN INDICES, not mm.

    Preserve original float32-MeV rounding and assignment. Duplicates, negative
    or fractional indices stop conversion instead of silently changing showers.
    Upper-outside hits are skipped as in the old converter and explicitly counted.
    """
    shower = np.zeros(shape, dtype=np.float32)
    seen = set()
    outside = 0
    outside_energy = 0.
    # Index the C++ collection explicitly. Never pass a PyROOT vector/proxy
    # to NumPy and do not construct a C++ iterator for each event.
    for hit_index in range(len(hits)):
        hit = hits[hit_index]
        position = hit.position
        xyz = (float(position.x), float(position.y), float(position.z))
        if not all(math.isfinite(a) for a in xyz):
            raise ValueError('Nonfinite voxel index')
        index = tuple(int(a) for a in xyz)
        if any(abs(a-b) > 1e-5 for a,b in zip(xyz,index)) or any(a<0 for a in index):
            raise ValueError('Negative/fractional voxel index: '+str(xyz))
        # Matches std::vector<float>(p.energy * 1.e3) in the original converter.
        energy = np.float32(float(hit.energy)*1000.)
        if not np.isfinite(energy) or energy < 0:
            raise ValueError('Negative/nonfinite/overflowing deposited energy')
        value = np.float32(float(energy)/sampling_fraction)
        if not np.isfinite(value):
            raise ValueError('Energy overflow after sampling-fraction division')
        if any(a >= b for a,b in zip(index,shape)):
            outside += 1
            outside_energy += float(value)
            continue
        if index in seen:
            raise ValueError('Duplicate voxel index '+str(index)+'; inspect before conversion')
        seen.add(index)
        shower[index] = value
    return shower, outside, outside_energy



def root_projector(ROOT, python_only=True):
    """Compatibility entry point; the scalar projector is now always selected.

    ROOT 6.28-era vector array interfaces may enter Cling on conversion. Removing
    that interface also avoids errors cascading from a failed buffer-pointer
    lookup. This does not disable PyROOT's own binding machinery.
    """
    return project_hits


def memory_snapshot():
    """Best-effort process diagnostics; no ROOT calls and no Linux requirement."""
    info = {}
    try:
        for line in Path('/proc/self/status').read_text().splitlines():
            key, _, value = line.partition(':')
            if key in ('VmPeak', 'VmSize', 'VmHWM', 'VmRSS'):
                info[key] = value.strip()
    except OSError:
        pass
    try:
        with open('/proc/self/maps') as f:
            info['mapped_regions'] = sum(1 for _ in f)
        info['vm_max_map_count'] = int(Path('/proc/sys/vm/max_map_count').read_text())
    except (OSError, ValueError):
        pass
    try:
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        info['address_space_soft_limit_bytes'] = None if soft == resource.RLIM_INFINITY else soft
        info['address_space_hard_limit_bytes'] = None if hard == resource.RLIM_INFINITY else hard
    except (ImportError, AttributeError, OSError, ValueError):
        pass
    return info


def log_memory(stage):
    info = memory_snapshot()
    print(stage + ': ' + json.dumps(info, sort_keys=True), flush=True)
    return info


def fingerprint(path):
    stat = path.stat()
    return dict(path=str(path.resolve()), size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def code_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def git_info():
    try:
        p = subprocess.run(['git','rev-parse','HEAD'], cwd=Path(__file__).parent,
                           capture_output=True, text=True, check=True)
        dirty = subprocess.run(['git','status','--porcelain'], cwd=Path(__file__).parent,
                               capture_output=True, text=True, check=True)
        return dict(commit=p.stdout.strip(), working_tree_dirty=bool(dirty.stdout.strip()))
    except (OSError,subprocess.CalledProcessError):
        return None


def make_datasets(h5, shape, batch):
    options = dict(maxshape=(None,), chunks=(batch,), dtype='f4',
                   compression='gzip', compression_opts=9)
    for name in ('incident_energy','incident_theta','incident_phi'):
        h5.create_dataset(name, shape=(0,), **options)
    h5.create_dataset('showers', shape=(0,)+shape, maxshape=(None,)+shape,
                      chunks=(batch,)+shape, dtype='f4', compression='gzip', compression_opts=9)


def append_batch(h5, grids, conditions, count):
    if not count:
        return
    start = h5['showers'].shape[0]
    stop = start+count
    for name in h5:
        h5[name].resize(stop, axis=0)
    h5['showers'][start:stop] = grids[:count]
    for col,name in enumerate(('incident_energy','incident_theta','incident_phi')):
        h5[name][start:stop] = conditions[:count,col]


def convert(args):
    source = Path(args.inputFile).expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError('Missing or empty ROOT file: '+str(source))
    outdir = Path(args.outputDir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    name = args.output_name or source.stem+'_v2.h5'
    if Path(name).name != name or not name.endswith('_v2.h5'):
        raise ValueError('--output-name must be a basename ending _v2.h5')
    out = outdir/name
    cache = out.with_suffix('.selection.npz')
    report = out.with_suffix('.json')
    before = fingerprint(source)
    g = GEOMETRIES[args.detector]
    settings = dict(detector=args.detector, geometry=g, geometry_source_commit=GEOMETRY_COMMIT,
                    shape=[args.numR,args.numPhi,args.numZ], collection=args.collection,
                    sampling_fraction=args.samplingFraction, hit_projection='python_scalar_indexed',
                    boundary_margin_mm=args.boundary_margin_mm, time_margin_ns=args.time_margin_ns,
                    script_sha256=code_hash())
    if any(p.exists() for p in (out,cache,report)):
        if args.skip_existing and all(p.is_file() for p in (out,cache,report)):
            old = json.loads(report.read_text())
            if old.get('status') == 'complete' and old.get('input') == before and old.get('settings') == settings:
                with h5py.File(out,'r') as f:
                    if bool(f.attrs.get('conversion_complete',False)) and f['showers'].shape[0]==old['new_kept']:
                        print('Already complete with identical input/settings: '+str(out),flush=True)
                        return
        raise FileExistsError('Output exists or changed configuration; use a NEW output directory: '+str(out))
    tmp = {p:p.with_name(p.name+'.partial') for p in (out,cache,report)}
    if any(p.exists() for p in tmp.values()):
        raise FileExistsError('A .partial file exists; inspect/remove failed temporary outputs before retrying')
    print('Converter ' + VERSION + ': scalar/indexed hit copying; no custom JIT helper', flush=True)
    log_memory('Before ROOT import')
    import ROOT
    ROOT.gROOT.SetBatch(True)
    if hasattr(ROOT,'IsImplicitMTEnabled') and ROOT.IsImplicitMTEnabled():
        ROOT.DisableImplicitMT()
    log_memory('After ROOT import')
    rf = ROOT.TFile.Open(str(source),'READ')
    if not rf or rf.IsZombie():
        raise OSError('ROOT could not open '+str(source))
    try:
        tree = rf.Get('events')
        if not tree or not tree.InheritsFrom('TTree'):
            raise ValueError('Expected a TTree named events')
        for branch in ('MCParticles',args.collection):
            if not tree.GetBranch(branch):
                raise ValueError('Required ROOT branch missing: '+branch)
        tree.SetBranchStatus('*',0)
        for branch in ('MCParticles',args.collection):
            tree.SetBranchStatus(branch,1)
            tree.SetBranchStatus(branch+'*',1)
        tree.SetCacheSize(8*1024*1024)
        projector = project_hits
        n = int(tree.GetEntries())
        if n == 0:
            raise ValueError('Empty input tree')
        d = {'root_row':np.arange(n,dtype=np.int64), 'nmc':np.zeros(n,dtype=np.int32),
             'reason':np.zeros(n,dtype=np.uint8), 'incident_energy_MeV':np.full(n,np.nan,dtype=np.float32),
             'incident_theta':np.full(n,np.nan,dtype=np.float32),
             'incident_phi':np.full(n,np.nan,dtype=np.float32),
             'etot_MeV':np.zeros(n,dtype=np.float64),
             'endpoint_x_mm':np.full(n,np.nan), 'endpoint_y_mm':np.full(n,np.nan),
             'endpoint_z_mm':np.full(n,np.nan)}
        shape = tuple(settings['shape'])
        grids = np.zeros((args.batch_size,)+shape,dtype=np.float32)
        conditions = np.zeros((args.batch_size,3),dtype=np.float32)
        used = 0; outside_count = 0; outside_energy = 0.; outside_events = 0
        with h5py.File(tmp[out],'x') as h5:
            make_datasets(h5,shape,args.batch_size)
            h5.attrs['conversion_complete'] = False
            h5.attrs['converter_version'] = VERSION
            h5.attrs['settings_json'] = json.dumps(settings,sort_keys=True)
            h5.attrs['source_root'] = str(source)
            h5.attrs['selection'] = 'upstream primary endpoint OR early stored secondary upstream'
            h5.attrs['units'] = 'energies: MeV; angles: rad'
            h5.attrs['scope'] = 'Original barrel scoring collection; approximate ECAL-entry selection'
            for row in range(n):
                if int(tree.GetEntry(row)) <= 0:
                    raise OSError('ROOT read failed at row '+str(row))
                try:
                    particles = getattr(tree,'MCParticles')
                    d['nmc'][row] = len(particles)
                    reason = classify(particles,g,args.boundary_margin_mm,args.time_margin_ns)
                    d['reason'][row] = reason
                    if len(particles):
                        p = particles[0]
                        px,py,pz = float(p.momentum.x),float(p.momentum.y),float(p.momentum.z)
                        norm = math.sqrt(px*px+py*py+pz*pz)
                        if norm>0 and math.isfinite(norm):
                            d['incident_energy_MeV'][row] = norm*1000.
                            d['incident_theta'][row] = math.atan2(math.hypot(px,py),pz)
                            d['incident_phi'][row] = math.atan2(py,px)
                        for axis in ('x','y','z'):
                            d['endpoint_'+axis+'_mm'][row] = float(getattr(p.endpoint,axis))
                    shower,outside,escaped = projector(getattr(tree,args.collection),shape,args.samplingFraction)
                    d['etot_MeV'][row] = shower.sum(dtype=np.float64)
                    outside_count += outside; outside_energy += escaped; outside_events += int(outside>0)
                    if reason == 0:
                        vals = [d[k][row] for k in ('incident_energy_MeV','incident_theta','incident_phi')]
                        if not all(np.isfinite(x) for x in vals) or vals[0] <= 0:
                            raise ValueError('Invalid stored conditioning at row '+str(row))
                        grids[used] = shower
                        conditions[used] = vals
                        used += 1
                        if used == args.batch_size:
                            append_batch(h5,grids,conditions,used); used = 0
                    if row == 0 or (row+1)%args.progress_every==0 or row+1==n:
                        print(f'{source.name}: {row+1}/{n}; written {h5["showers"].shape[0]+used}',flush=True)
                        log_memory('Process resources')
                    del particles, shower
                except Exception as exc:
                    raise RuntimeError(
                        f'Conversion failed at ROOT row {row} (0-based), source {source.name}: {exc}'
                    ) from exc
            append_batch(h5,grids,conditions,used)
            new = d['reason']==0; old = d['nmc']==1
            if h5['showers'].shape[0] != int(new.sum()):
                raise RuntimeError('Internal count mismatch')
            if fingerprint(source) != before:
                raise RuntimeError('ROOT input changed during conversion')
            summary = dict(status='complete', version=VERSION,
                created_utc=datetime.now(timezone.utc).isoformat(), input=before, settings=settings,
                git=git_info(), resources_at_end=memory_snapshot(), software=dict(ROOT=ROOT.gROOT.GetVersion(),numpy=np.__version__,
                h5py=h5py.__version__,python=sys.version.split()[0]),
                output_file=str(out), diagnostic_cache=str(cache), total_events=n,
                new_kept=int(new.sum()), old_numMC1_kept=int(old.sum()),
                counts_by_reason={REASONS[k]:int((d['reason']==k).sum()) for k in REASONS},
                overlap=dict(both_keep=int((new&old).sum()),new_kept_old_rejected=int((new&~old).sum()),
                             new_rejected_old_kept=int((~new&old).sum()),both_rejected=int((~new&~old).sum())),
                upper_outside_hits_skipped=outside_count, upper_outside_events=outside_events,
                upper_outside_energy_MeV=outside_energy,
                zero_energy_kept=int((new&(d['etot_MeV']==0)).sum()),
                limitation='Geometry/nominal-time truth approximation, NOT measured entry state; may miss scattering without saved recoil or misclassify returning primaries.')
            h5.attrs['report_json'] = json.dumps(summary,sort_keys=True)
            h5.attrs['conversion_complete'] = True
        # Read back small headers/incident arrays, not all shower tensors.
        with h5py.File(tmp[out],'r') as h5:
            for key in ('showers','incident_energy','incident_theta','incident_phi'):
                ds=h5[key]
                if ds.shape[0]!=summary['new_kept'] or ds.dtype!=np.dtype('float32'):
                    raise RuntimeError('Output schema check failed: '+key)
        with open(tmp[cache],'xb') as f:
            np.savez_compressed(f,**d,metadata_json=np.array(json.dumps(summary)))
        with open(tmp[report],'x') as f:
            json.dump(summary,f,indent=2); f.write('\n')
        # Publish the HDF5 last; .partial is never a finished dataset.
        os.replace(tmp[cache],cache); os.replace(tmp[report],report); os.replace(tmp[out],out)
        print(json.dumps({k:summary[k] for k in ('output_file','total_events','new_kept','old_numMC1_kept','counts_by_reason','overlap')},indent=2))
        if outside_count:
            print('REVIEW: upper-outside hits were skipped as in the original converter; see JSON',file=sys.stderr)
    finally:
        rf.Close()


def parse_args(argv=None):
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i','--inputFile',required=True)
    p.add_argument('-o','--outputDir',required=True)
    p.add_argument('--detector',choices=sorted(GEOMETRIES),required=True)
    p.add_argument('--output-name',help='e.g. part1_v2.h5; default original_basename_v2.h5')
    p.add_argument('--collection',default='ECalBarrelCollection')
    p.add_argument('--numR',type=int,default=9); p.add_argument('--numPhi',type=int,default=16)
    p.add_argument('--numZ',type=int,default=45)
    p.add_argument('--samplingFraction',type=float,default=1.)
    p.add_argument('--boundary-margin-mm',type=float,default=1.)
    p.add_argument('--time-margin-ns',type=float,default=.1)
    p.add_argument('--batch-size',type=int,default=64)
    p.add_argument('--python-hits',action='store_true',default=True,
                   help='Compatibility flag: scalar/indexed hit copying is now always used')
    p.add_argument('--progress-every',type=int,default=1000,
                   help='Print progress and process memory diagnostics every N events (default 1000)')
    p.add_argument('--skip-existing',action='store_true',help='Skip ONLY a complete output with identical input and code/settings')
    a=p.parse_args(argv)
    if min(a.numR,a.numPhi,a.numZ,a.batch_size,a.progress_every)<=0: p.error('Dimensions/batch/progress interval must be positive')
    if not math.isfinite(a.samplingFraction) or a.samplingFraction<=0: p.error('Invalid samplingFraction')
    if not 0<=a.boundary_margin_mm<100 or not 0<=a.time_margin_ns<100: p.error('Invalid margins')
    return a


if __name__=='__main__':
    try:
        convert(parse_args())
    except Exception as e:
        print('ERROR: '+str(e),file=sys.stderr)
        print('Process resources: '+json.dumps(memory_snapshot(),sort_keys=True),file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

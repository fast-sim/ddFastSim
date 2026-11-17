"""No ROOT installation required: math, voxel conversion, HDF5, and mock-reader tests.

These are not a substitute for a test on the actual EDM4hep ROOT production.
"""
import contextlib
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import subprocess
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch

import h5py
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1] / 'scripts'))
import edm4hep2h5_v2 as converter


def vec(x=0.,y=0.,z=0.): return NS(x=x,y=y,z=z)
def photon(endpoint=(2200.,0.,0.),momentum=(50.,0.,0.)):
    return NS(PDG=22,generatorStatus=1,vertex=vec(),endpoint=vec(*endpoint),
              momentum=vec(*momentum),time=0.)
def secondary(vertex=(100.,0.,0.),time=1.):
    return NS(PDG=11,generatorStatus=0,vertex=vec(*vertex),endpoint=vec(*vertex),
              momentum=vec(.001),time=time)
def hit(index=(0,0,0),energy=.001): return NS(position=vec(*index),energy=energy)


class FakeTree:
    def __init__(self,rows): self.rows=rows
    def InheritsFrom(self,name): return name=='TTree'
    def GetEntries(self): return len(self.rows)
    def GetBranch(self,name): return name in ('MCParticles','ECalBarrelCollection')
    def SetBranchStatus(self,*args): pass
    def SetCacheSize(self,*args): pass
    def GetEntry(self,i):
        if i>=len(self.rows): return 0
        self.MCParticles,self.ECalBarrelCollection=self.rows[i]
        return 1


class FakeFile:
    def __init__(self,rows): self.tree=FakeTree(rows)
    def IsZombie(self): return False
    def Get(self,name): return self.tree if name=='events' else None
    def Close(self): pass


def fake_root(rows):
    return NS(gROOT=NS(SetBatch=lambda *_:None,GetVersion=lambda:'MOCK-NOT-REAL-ROOT'),
              IsImplicitMTEnabled=lambda:False,TFile=NS(Open=lambda *_:FakeFile(rows)))


class GeometryTests(unittest.TestCase):
    def test_CLD_face_and_corner(self):
        g=converter.GEOMETRIES['CLD']
        self.assertAlmostEqual(converter.support(2150,0,g),2150)
        phi=math.pi/12
        r=2150/math.cos(phi)
        self.assertAlmostEqual(converter.support(r*math.cos(phi),r*math.sin(phi),g),2150)
        self.assertAlmostEqual(converter.entry_distance(math.cos(phi),math.sin(phi),0,g),r)

    def test_cld_primary_endpoint(self):
        g=converter.GEOMETRIES['CLD']
        self.assertEqual(converter.classify([photon((1000,0,0))],g),3)
        self.assertEqual(converter.classify([photon((2150,0,0))],g),0)
        self.assertEqual(converter.classify([photon((2149.5,0,0))],g),0)

    def test_secondary_time_and_location(self):
        g=converter.GEOMETRIES['CLD']
        self.assertEqual(converter.classify([photon(),secondary(time=1)],g),4)
        self.assertEqual(converter.classify([photon(),secondary(time=400)],g),0)
        self.assertEqual(converter.classify([photon(),secondary(vertex=(2200,0,0),time=1)],g),0)
        # Equality at time margin is deliberately not rejected.
        t=2150/converter.C-.1
        self.assertEqual(converter.classify([photon(),secondary(time=t)],g),0)

    def test_allegro_region_not_active_layer(self):
        g=converter.GEOMETRIES['ALLEGRO']
        self.assertAlmostEqual(converter.entry_distance(1,0,0,g),2149)
        self.assertEqual(converter.classify([photon((2140,0,0))],g),3)
        self.assertEqual(converter.classify([photon((2155,0,0))],g),0)
        self.assertEqual(g['endcap_z'],3260)
        self.assertEqual(g['endcap_inner_r'],410)

    def test_endcap_face_and_hole(self):
        for g in converter.GEOMETRIES.values():
            self.assertTrue(math.isinf(converter.entry_distance(0,0,1,g)))
            px=1000/g['endcap_z']
            expected=math.hypot(1000,g['endcap_z'])
            self.assertAlmostEqual(converter.entry_distance(px,0,1,g),expected)

    def test_primary_identity(self):
        g=converter.GEOMETRIES['CLD']
        self.assertEqual(converter.classify([],g),1)
        p=photon();p.generatorStatus=0
        self.assertEqual(converter.classify([p],g),1)
        p=photon();p.PDG=11
        self.assertEqual(converter.classify([p],g),1)
        p=photon();p.vertex.x=2.
        self.assertEqual(converter.classify([p],g),1)
        self.assertEqual(converter.classify([photon(),photon()],g),1)

    def test_nonfinite_not_accepted(self):
        p=photon();p.endpoint.x=math.nan
        with self.assertRaises(ValueError): converter.classify([p],converter.GEOMETRIES['CLD'])


class ProjectionTests(unittest.TestCase):
    def test_units_order_final_bin_zero(self):
        a,n,e=converter.project_hits([hit((8,15,44),.002),hit((0,0,1),.001)],(9,16,45),1.)
        self.assertEqual(a[8,15,44],2.);self.assertEqual(a[0,0,1],1.)
        self.assertEqual(a.sum(),3.);self.assertEqual(a.dtype,np.float32)
        self.assertEqual(n,0)
        a,_,_=converter.project_hits([],(9,16,45),1.)
        self.assertEqual(a.sum(),0)

    def test_original_rounding_and_scaling(self):
        x=.0123456789
        a,_,_=converter.project_hits([hit(energy=x)],(9,16,45),.031)
        expected=np.float32(float(np.float32(x*1000.))/.031)
        self.assertEqual(a[0,0,0],expected)

    def test_upper_outside_accounted(self):
        a,n,e=converter.project_hits([hit((9,0,0),.003)],(9,16,45),1.)
        self.assertEqual((n,e,a.sum()),(1,3.,0.))

    def test_bad_indices_duplicates_or_energy_stop(self):
        for hits in ([hit((-1,0,0))],[hit((.5,0,0))],[hit(),hit()],
                     [hit(energy=-1)],[hit(energy=float('nan'))]):
            with self.assertRaises(ValueError): converter.project_hits(hits,(9,16,45),1.)


class ConversionTests(unittest.TestCase):
    def test_stream_order_readback_cache_and_safe_restart(self):
        rows=[([photon()],[hit((8,15,44),.002)]),
              ([photon((1000,0,0))],[hit(energy=.003)]),
              ([photon(),secondary(time=400)],[hit(energy=.004)]),
              ([photon(),secondary(time=1)],[hit(energy=.005)]),
              ([photon()],[])]
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)/'source.root';root.write_bytes(b'mock input')
            out=Path(temp)/'out'
            a=converter.parse_args(['-i',str(root),'-o',str(out),'--detector','CLD',
                                   '--output-name','part1_v2.h5','--batch-size','2','--python-hits'])
            before=root.read_bytes()
            with patch.dict(sys.modules,{'ROOT':fake_root(rows)}), contextlib.redirect_stdout(io.StringIO()):
                converter.convert(a)
            self.assertEqual(root.read_bytes(),before)
            with h5py.File(out/'part1_v2.h5','r') as f:
                self.assertEqual(set(f),{'showers','incident_energy','incident_theta','incident_phi'})
                self.assertEqual(f['showers'].shape,(3,9,16,45))
                np.testing.assert_array_equal(f['showers'][:].sum((1,2,3)),[2.,4.,0.])
                self.assertEqual(f['showers'].compression,'gzip')
                self.assertEqual(f['showers'].compression_opts,9)
                self.assertTrue(f.attrs['conversion_complete'])
            with np.load(out/'part1_v2.selection.npz',allow_pickle=False) as f:
                np.testing.assert_array_equal(f['reason'],[0,3,0,4,0])
                np.testing.assert_array_equal(f['etot_MeV'],[2,3,4,5,0])
                np.testing.assert_array_equal(f['root_row'],np.arange(5))
            report=json.loads((out/'part1_v2.json').read_text())
            self.assertEqual(report['overlap'],{'both_keep':2,'new_kept_old_rejected':1,'new_rejected_old_kept':1,'both_rejected':1})
            with self.assertRaises(FileExistsError): converter.convert(a)
            a.skip_existing=True
            with contextlib.redirect_stdout(io.StringIO()): converter.convert(a)
            a.time_margin_ns=.2
            with self.assertRaises(FileExistsError): converter.convert(a)

    def test_fail_leaves_partial_not_finished_output(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'source.root';p.write_bytes(b'mock')
            a=converter.parse_args(['-i',str(p),'-o',str(Path(t)/'out'),
                                   '--detector','ALLEGRO','--output-name','part1_v2.h5','--python-hits'])
            with patch.dict(sys.modules,{'ROOT':fake_root([([photon()],[hit(),hit()])])}):
                with self.assertRaisesRegex(RuntimeError, "ROOT row 0.*Duplicate voxel"):
                    converter.convert(a)
            self.assertFalse((Path(t)/'out/part1_v2.h5').exists())
            self.assertTrue((Path(t)/'out/part1_v2.h5.partial').exists())


class RepositoryLayoutTests(unittest.TestCase):
    def test_plot_cli_from_unrelated_working_directory(self):
        plot = Path(__file__).resolve().parents[1] / 'zenodo_utils' / 'plot_selection_v2.py'
        with tempfile.TemporaryDirectory() as work:
            result = subprocess.run([sys.executable, str(plot), '--help'], cwd=work,
                                    capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--outdir', result.stdout)



class ScalarOnlySequence:
    """A proxy that fails on iteration or array-interface access."""
    def __init__(self, values): self.values = values
    def __len__(self): return len(self.values)
    def __getitem__(self, i): return self.values[i]
    def __iter__(self): raise AssertionError('Must use indexed scalar reads')
    @property
    def __array_interface__(self):
        raise AssertionError('Must not request a C++/NumPy array interface')


class StreamingRegressionTests(unittest.TestCase):
    def test_default_never_uses_interpreter(self):
        class ForbiddenROOT:
            def __getattr__(self, name):
                raise AssertionError('ROOT helper unexpectedly used: ' + name)
        self.assertIs(converter.root_projector(ForbiddenROOT(), False), converter.project_hits)

    def test_scalar_proxy_only(self):
        hits=ScalarOnlySequence([hit((8,15,44),.002),hit((0,0,0),.001)])
        grid,n,e=converter.project_hits(hits,(9,16,45),1.)
        self.assertEqual((grid.sum(),n,e),(3.,0,0.))
        particles=ScalarOnlySequence([photon(),secondary(time=400)])
        self.assertEqual(converter.classify(particles,converter.GEOMETRIES['CLD']),0)

    def test_cli_default_and_compatibility_flag(self):
        for extra in ([], ['--python-hits']):
            args=converter.parse_args(['-i','x.root','-o','out','--detector','CLD']+extra)
            self.assertTrue(args.python_hits)
            self.assertEqual(args.batch_size,64)

    def test_conversion_works_without_opt_in_flag_or_interpreter(self):
        rows=[(ScalarOnlySequence([photon()]),ScalarOnlySequence([hit((8,15,44),.002)]))]
        with tempfile.TemporaryDirectory() as t:
            src=Path(t)/'source.root';src.write_bytes(b'mock')
            args=converter.parse_args(['-i',str(src),'-o',str(Path(t)/'out'),
                                      '--detector','CLD','--output-name','part1_v2.h5'])
            with patch.dict(sys.modules,{'ROOT':fake_root(rows)}), contextlib.redirect_stdout(io.StringIO()):
                converter.convert(args)
            with h5py.File(Path(t)/'out/part1_v2.h5','r') as f:
                self.assertEqual(f['showers'][0,8,15,44],2.)
                self.assertTrue(f.attrs['conversion_complete'])
            report=json.loads((Path(t)/'out/part1_v2.json').read_text())
            self.assertEqual(report['version'],'lemurs-v2-converter-1.0.1')
            self.assertEqual(report['settings']['hit_projection'],'python_scalar_indexed')

    def test_resource_diagnostics_serializable(self):
        json.dumps(converter.memory_snapshot(),allow_nan=False)

if __name__=='__main__': unittest.main()

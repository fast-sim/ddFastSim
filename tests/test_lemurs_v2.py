"""Tests for the RDataFrame converter; no ROOT installation is required."""

import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import edm4hep2h5_v2 as converter


def object_column(*rows):
    result = np.empty(len(rows), dtype=object)
    for index, row in enumerate(rows):
        result[index] = np.asarray(row)
    return result


def sample_arrays():
    return {
        "energyCell": object_column([2., 1.], [4.], []),
        "rhoCell": object_column([8, 0], [9], []),
        "phiCell": object_column([15, 0], [0], []),
        "zCell": object_column([44, 1], [0], []),
        "energyMC": np.asarray([50000., 60000., 70000.]),
        "thetaMC": np.asarray([1., 1.1, 1.2]),
        "phiMC": np.asarray([0., .1, .2]),
    }


class FakeInterpreter:
    def __init__(self):
        self.code = None

    def Declare(self, code):
        self.code = code
        return True


class FakeDataFrame:
    def __init__(self, arrays):
        self.arrays = arrays
        self.filter = None
        self.defines = []

    def Filter(self, expression, label):
        self.filter = (expression, label)
        return self

    def Define(self, name, expression):
        self.defines.append((name, expression))
        return self

    def AsNumpy(self, columns):
        if set(columns) != set(converter.COLUMNS):
            raise AssertionError(columns)
        return self.arrays


class FakeROOT:
    def __init__(self, arrays):
        self.gInterpreter = FakeInterpreter()
        self.frame = FakeDataFrame(arrays)

    def RDataFrame(self, tree, source):
        if tree != "events":
            raise AssertionError(tree)
        return self.frame

    @staticmethod
    def IsImplicitMTEnabled():
        return False


class GeometryTests(unittest.TestCase):
    def test_support_coordinates(self):
        cld = converter.GEOMETRIES["CLD"]
        allegro = converter.GEOMETRIES["ALLEGRO"]
        self.assertAlmostEqual(converter.support(2150., 0., cld), 2150.)
        self.assertAlmostEqual(converter.support(3., 4., allegro), 5.)

    def test_filter_uses_detector_geometry(self):
        args = converter.parse_args([
            "-i", "input.root", "-o", "output", "--detector", "CLD"])
        expression = converter.selection_expression(args)
        self.assertIn("keep_event(MCParticles", expression)
        self.assertIn("2150", expression)
        self.assertIn("2307", expression)


class DataFrameTests(unittest.TestCase):
    def test_root_style_filter_define_chain(self):
        root = FakeROOT(sample_arrays())
        args = converter.parse_args([
            "-i", "input.root", "-o", "output", "--detector", "ALLEGRO"])
        frame = converter.build_dataframe(root, Path("input.root"), args)
        self.assertIs(frame, root.frame)
        self.assertIn("keep_event(MCParticles", frame.filter[0])
        self.assertEqual(frame.filter[1], "pre-ECAL interaction veto")
        self.assertEqual([name for name, _ in frame.defines],
                         list(converter.COLUMNS))
        self.assertIn("bool keep_event", root.gInterpreter.code)

    def test_shower_batch_scaling_and_outside_hit(self):
        block, skipped = converter.shower_block(
            sample_arrays(), 0, 3, (9, 16, 45), 2.)
        self.assertEqual(block.shape, (3, 9, 16, 45))
        self.assertEqual(block[0, 8, 15, 44], 1.)
        self.assertEqual(block[0, 0, 0, 1], .5)
        self.assertEqual(block[1].sum(), 0.)
        self.assertEqual(skipped, 1)


class ConversionTests(unittest.TestCase):
    def test_hdf5_schema_and_atomic_publish(self):
        root = FakeROOT(sample_arrays())
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source.root"
            source.write_bytes(b"mock ROOT")
            output_dir = Path(temporary) / "output"
            args = converter.parse_args([
                "-i", str(source), "-o", str(output_dir), "--detector", "CLD",
                "--output-name", "part1_v2.h5", "--batch-size", "2"])
            with patch.dict(sys.modules, {"ROOT": root}), \
                    contextlib.redirect_stdout(io.StringIO()):
                output = converter.convert(args)

            self.assertEqual(output, output_dir / "part1_v2.h5")
            self.assertFalse((output_dir / "part1_v2.h5.partial").exists())
            self.assertFalse((output_dir / "part1_v2.json").exists())
            self.assertFalse((output_dir / "part1_v2.selection.npz").exists())
            with h5py.File(output, "r") as h5:
                self.assertEqual(set(h5), {
                    "showers", "incident_energy", "incident_theta",
                    "incident_phi"})
                self.assertEqual(h5["showers"].shape, (3, 9, 16, 45))
                self.assertEqual(h5["showers"].dtype, np.dtype("float32"))
                self.assertEqual(h5["showers"].compression, "gzip")
                self.assertEqual(h5["showers"].compression_opts, 9)
                self.assertEqual(h5["showers"][0].sum(), 3.)
                self.assertEqual(h5["showers"][1].sum(), 0.)
                self.assertEqual(h5.attrs["detector"], "CLD")

            args.skip_existing = True
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(converter.convert(args), output)


class RepositoryLayoutTests(unittest.TestCase):
    def test_plot_cli_from_unrelated_working_directory(self):
        plot = (Path(__file__).resolve().parents[1] /
                "zenodo_utils" / "plot_selection_v2.py")
        with tempfile.TemporaryDirectory() as work:
            result = subprocess.run(
                [sys.executable, str(plot), "--help"], cwd=work,
                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--outdir", result.stdout)


if __name__ == "__main__":
    unittest.main()

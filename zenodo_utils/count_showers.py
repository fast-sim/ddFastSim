# count_showers.py
from pathlib import Path
import h5py

DETECTORS = ["Par04SiW", "Par04SciPb", "ODD", "FCCeeCLD", "FCCeeALLEGRO"]

def try_len_of_dataset(h5, name):
    if name in h5:
        d = h5[name]
        try:
            return len(d)
        except TypeError:
            pass
    # search case-insensitively
    for k in h5.keys():
        if name.lower() in k.lower():
            try:
                return len(h5[k])
            except TypeError:
                continue
    return None

def generic_count(h5):
    # As a fallback: pick the first dataset with a sensible leading dimension.
    def walk(g):
        for k, v in g.items():
            if isinstance(v, h5py.Dataset):
                shape = v.shape
                if shape and shape[0] > 0:
                    return len(v)
            elif isinstance(v, h5py.Group):
                out = walk(v)
                if out is not None:
                    return out
        return None
    return walk(h5)

def count_showers(path: Path):
    try:
        with h5py.File(path, "r") as h5:
            # 1) explicit dataset 'showers'
            n = try_len_of_dataset(h5, "showers")
            if n is not None:
                return n
            # 2) attribute 'nShowers'
            for node in [h5] + [h5[k] for k in h5.keys()]:
                if hasattr(node, "attrs") and "nShowers" in node.attrs:
                    return int(node.attrs["nShowers"])
            # 3) common event datasets
            for candidate in ("Events", "events", "event_data"):
                n = try_len_of_dataset(h5, candidate)
                if n is not None:
                    return n
            # 4) fallback heuristic
            n = generic_count(h5)
            return n
    except Exception as e:
        return f"ERROR: {e}"

def main(root="to_upload/"):
    root = Path(root).resolve()
    for det in DETECTORS:
        for split in ("training", "testing"):
            folder = root / det / split
            if not folder.exists():
                continue
            for f in sorted(folder.glob("*.h5")):
                n = count_showers(f)
                print(f"{f.name}: {n}")

if __name__ == "__main__":
    main()

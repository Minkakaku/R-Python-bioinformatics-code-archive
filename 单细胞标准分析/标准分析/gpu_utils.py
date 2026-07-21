import time
import warnings


def safe_import_rapids_singlecell():
    try:
        import rapids_singlecell as rsc

        return rsc
    except Exception as exc:
        warnings.warn(f"rapids_singlecell import failed: {exc}")
        return None


def check_gpu_available(verbose=True):
    try:
        import cupy as cp
        import rapids_singlecell  # noqa: F401

        n_gpu = cp.cuda.runtime.getDeviceCount()
        if n_gpu < 1:
            if verbose:
                print("[GPU] unavailable, no CUDA device found.")
            return False
        name = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
        if verbose:
            print(f"[GPU] available: {name}")
        return True
    except Exception as exc:
        if verbose:
            print(f"[GPU] unavailable, fallback to CPU. Reason: {exc}")
        return False


def to_gpu(adata):
    rsc = safe_import_rapids_singlecell()
    if rsc is None:
        raise RuntimeError("rapids_singlecell is not available")
    rsc.get.anndata_to_GPU(adata)
    return adata


def to_cpu(adata):
    rsc = safe_import_rapids_singlecell()
    if rsc is None:
        return adata
    try:
        rsc.get.anndata_to_CPU(adata)
    except Exception:
        pass
    return adata


class Timer:
    def __init__(self, name):
        self.name = name
        self.start = None

    def __enter__(self):
        self.start = time.time()
        print(f"[START] {self.name}")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        elapsed = time.time() - self.start
        print(f"[DONE] {self.name}: {elapsed:.2f} sec")

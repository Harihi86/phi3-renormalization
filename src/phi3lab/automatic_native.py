"""Optional host-checked generic C++ multiprecision BPH executor.

The ordinary Python compiler remains sufficient without this module's optional
library. Build explicitly with tools/build_automatic_native.py and external
Boost headers. No download or compilation occurs when importing this module.
"""
import ctypes
from collections import OrderedDict
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import weakref
from functools import lru_cache

from .automatic_bph import _external_gram


class NativeBackendUnavailable(RuntimeError):
    """The optional library is missing, stale, or built for another host."""


_ROOT = Path(__file__).resolve().parents[1]
_BUILD = _ROOT / 'native' / 'build'
_EXECUTORS = OrderedDict()
_LOADED_LIBRARY_SIGNATURE = None


def _library_path():
    return _BUILD / ('automatic_bph.dll' if sys.platform == 'win32' else
                     'libautomatic_bph.dylib' if sys.platform == 'darwin' else
                     'libautomatic_bph.so')


@lru_cache(maxsize=32)
def _hash_at_stat(path, stat_key):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _file_hash(path):
    stat=path.stat()
    return _hash_at_stat(str(path),(stat.st_dev,stat.st_ino,stat.st_size,
                                   stat.st_mtime_ns,stat.st_ctime_ns))


def _verified_metadata():
    if os.environ.get('PHI3LAB_DISABLE_NATIVE') == '1':
        raise NativeBackendUnavailable('Optional BPH executor disabled by PHI3LAB_DISABLE_NATIVE=1')
    path = _library_path()
    try:
        metadata = json.loads((_BUILD / 'automatic_native_build.json').read_text())
        if metadata['host'] != {'system': platform.system(), 'machine': platform.machine()}:
            raise NativeBackendUnavailable('Optional BPH library belongs to another host; rebuild it')
        source = _ROOT / 'native' / 'automatic_bph.cpp'
        if metadata['source_sha256'] != _file_hash(source):
            raise NativeBackendUnavailable('Optional BPH library is stale; rebuild it')
        if metadata['library_sha256'] != _file_hash(path):
            raise NativeBackendUnavailable('Optional BPH library does not match its build record')
    except NativeBackendUnavailable:
        raise
    except (OSError, ValueError, KeyError) as exc:
        raise NativeBackendUnavailable('Build the optional multiprecision library with tools/build_automatic_native.py') from exc
    return metadata


def _load_library():
    global _LOADED_LIBRARY_SIGNATURE
    metadata=_verified_metadata()
    signature=(metadata['source_sha256'],metadata['library_sha256'])
    if _LOADED_LIBRARY_SIGNATURE is not None and _LOADED_LIBRARY_SIGNATURE != signature:
        raise NativeBackendUnavailable('Optional BPH library changed after loading; restart this process')
    try:
        library=ctypes.CDLL(str(_library_path()))
    except OSError as exc:
        raise NativeBackendUnavailable('Optional BPH library could not be loaded on this host') from exc
    pointer = ctypes.c_void_p
    doubles = ctypes.POINTER(ctypes.c_double)
    integers = ctypes.POINTER(ctypes.c_int)
    library.phi3lab_bph_abi_version.restype = ctypes.c_int
    if library.phi3lab_bph_abi_version() != 1:
        raise NativeBackendUnavailable('Optional BPH library ABI mismatch')
    _LOADED_LIBRARY_SIGNATURE=signature
    library.phi3lab_bph_last_error.restype = ctypes.c_char_p
    library.phi3lab_bph_create.argtypes = (integers, ctypes.c_size_t)
    library.phi3lab_bph_create.restype = pointer
    library.phi3lab_bph_destroy.argtypes = (pointer,)
    library.phi3lab_bph_batch.argtypes = (pointer, ctypes.c_size_t, doubles, doubles,
                                        ctypes.c_double, doubles, integers, doubles)
    library.phi3lab_bph_batch.restype = ctypes.c_int
    library.phi3lab_bph_batch_parallel.argtypes = (pointer, ctypes.c_size_t, doubles,
                                                  doubles, ctypes.c_double, ctypes.c_int,
                                                  doubles, integers, doubles)
    library.phi3lab_bph_batch_parallel.restype = ctypes.c_int
    library.phi3lab_bph_text.argtypes = (pointer, doubles, doubles, ctypes.c_double,
                                       ctypes.c_int, ctypes.c_char_p, ctypes.c_size_t)
    library.phi3lab_bph_text.restype = ctypes.c_int
    return library, metadata


def native_available():
    """Read-only availability check; the absence of Boost is harmless at run time."""
    try:
        _load_library()
        return True
    except NativeBackendUnavailable:
        return False


def _flatten_program(program):
    if program.get('schema') != 'phi3lab.recursive-gaussian-bph.v1':
        raise ValueError('Expected a compiled recursive Gaussian BPH program')
    loops, edges, legs = (program[k] for k in ('loops', 'edges', 'external_legs'))
    width = loops + legs - 1
    terms = program['terms']
    data = [1, loops, edges, legs, len(terms)]
    for term in terms:
        orders = term['jet_orders']
        data.extend((term['sign'], len(orders)))
        data.extend(orders)
        rows = term['scaled_routing']
        if len(rows) != edges or any(len(row) != width for row in rows):
            raise ValueError('Compiled routing shape does not match graph dimensions')
        for row in rows:
            for entry in row:
                if len(entry['powers']) != len(orders):
                    raise ValueError('Compiled routing monomial has incorrect dimension')
                data.append(entry['coefficient'])
                data.extend(entry['powers'])
    if any(type(value) is not int or not -(2**31) <= value < 2**31 for value in data):
        raise ValueError('Native program fields must be bounded integers')
    return data


class NativeBPHExecutor:
    """Execute CompiledBPH.to_dict(), with 50/100-digit checking and fallback.

    A library handle is immutable, so independent calls are safe. Batch outputs
    are Python lists and require no NumPy dependency. Pass an actual CompiledBPH
    instance to permit its adaptive Decimal fallback; serialized programs can be
    evaluated with fallback=False or raise if native precision is insufficient.
    """
    def __init__(self, compiled):
        self.compiled = compiled if hasattr(compiled, 'to_dict') else None
        self.program = compiled.to_dict() if self.compiled is not None else compiled
        data = _flatten_program(self.program)
        self.library, self.build_metadata = _load_library()
        values = (ctypes.c_int * len(data))(*data)
        self.handle = self.library.phi3lab_bph_create(values, len(data))
        if not self.handle:
            raise ValueError(self.library.phi3lab_bph_last_error().decode())
        self._finalizer = weakref.finalize(self, self.library.phi3lab_bph_destroy, self.handle)
        self.program_sha256 = hashlib.sha256(json.dumps(self.program, sort_keys=True,
                                                       separators=(',', ':')).encode()).hexdigest()

    def close(self):
        self._finalizer()
        self.handle = None

    def __enter__(self):
        return self

    def __exit__(self, *unused):
        self.close()

    def _inputs(self, alpha, invariants, m2):
        if not self.handle:
            raise RuntimeError('Native executor is closed')
        rows = [tuple(float(a) for a in row) for row in alpha]
        edges = self.program['edges']
        if any(len(row) != edges or any(not math.isfinite(a) or a <= 0 for a in row)
               or abs(math.fsum(row) - 1) > 1e-10 for row in rows):
            raise ValueError('One positive projective parameter per edge is required')
        if not math.isfinite(float(m2)) or m2 <= 0:
            raise ValueError('m2 must be positive and finite')
        if self.compiled is not None:
            graph = self.compiled.graph
        else:
            from .graphs import Graph
            graph = Graph.from_dict(self.program['graph'])
        gram = _external_gram(graph, invariants)
        flattened = [a for row in rows for a in row]
        aa = (ctypes.c_double * len(flattened))(*flattened)
        gg = (ctypes.c_double * max(1, len(gram)**2))(*(a for row in gram for a in row))
        return rows, aa, gg

    def batch_density(self, alpha, invariants=None, m2=1, *, diagnostics=None,
                      fallback=True, threads=None):
        rows, aa, gram = self._inputs(alpha, invariants, m2)
        count = len(rows)
        output = (ctypes.c_double * count)()
        statuses = (ctypes.c_int * count)()
        checks = (ctypes.c_double * count)()
        if threads is None:
            threads=int(os.environ.get('PHI3LAB_NATIVE_THREADS',min(4,os.cpu_count() or 1)))
        if type(threads) is not int or not 1<=threads<=64:
            raise ValueError('Use1--64 native threads')
        used_threads=min(threads,count) if count>=64 else 1
        code = self.library.phi3lab_bph_batch_parallel(self.handle, count, aa, gram, float(m2),
                                                     used_threads,output, statuses, checks)
        if code:
            raise ArithmeticError(self.library.phi3lab_bph_last_error().decode())
        result = list(output)
        failed = [i for i, status in enumerate(statuses) if status != 0]
        digits = 0
        if failed and fallback:
            if self.compiled is None:
                raise ArithmeticError('Native precision failed; a CompiledBPH instance is needed for Decimal fallback')
            for i in failed:
                value, used = self.compiled._adaptive_decimal(rows[i], invariants, m2)
                result[i] = float(value)
                digits = max(digits, used)
        elif failed:
            # Never return an unchecked value as an accepted result.
            for i in failed:
                result[i] = math.nan
        if diagnostics is not None:
            diagnostics.update(native_checked_points=count, native_precision_pair=[50, 100],
                               native_threads=used_threads,
                               native_status=list(statuses), native_relative_checks=list(checks),
                               native_precision_failed_points=len(failed),
                               native_double_underflow_points=sum(status == 3 for status in statuses),
                               decimal_fallback_points=len(failed) if fallback else 0,
                               maximum_decimal_digits=digits,
                               native_program_sha256=self.program_sha256,
                               native_library_sha256=self.build_metadata['library_sha256'],
                               numerical_precision_check_is_rigorous_error_bound=False)
        return result

    def density_text(self, alpha, invariants=None, m2=1, *, digits=100):
        """Raw 50/100-digit output for independent Decimal/mpmath verification."""
        if digits not in (50, 100):
            raise ValueError('Native text evaluation supports 50 or 100 digits')
        _, aa, gram = self._inputs([alpha], invariants, m2)
        output = ctypes.create_string_buffer(4096)
        if self.library.phi3lab_bph_text(self.handle, aa, gram, float(m2), digits,
                                        output, len(output)):
            raise ArithmeticError(self.library.phi3lab_bph_last_error().decode())
        return output.value.decode('ascii')


def try_native_batch_density(compiled, alpha, invariants=None, m2=1, *, diagnostics=None):
    """Native batch or None if unavailable; safe inside the Python evaluator.

    Never calls public density/batch_density on availability failure. Individual
    native precision failures still go directly to adaptive Decimal. This prevents
    recursion when the standard compiler invokes this entry for flagged rows.
    """
    key = hashlib.sha256(json.dumps(compiled.to_dict(), sort_keys=True,
                                    separators=(',', ':')).encode()).hexdigest()
    try:
        if os.environ.get('PHI3LAB_DISABLE_NATIVE') == '1':
            raise NativeBackendUnavailable('Optional BPH executor disabled by PHI3LAB_DISABLE_NATIVE=1')
        executor = _EXECUTORS.get(key)
        if executor is None:
            executor = NativeBPHExecutor(compiled)
            _EXECUTORS[key] = executor
            if len(_EXECUTORS) > 16:
                _EXECUTORS.popitem(last=False)
        else:
            current=_verified_metadata()
            if (current['source_sha256'] != executor.build_metadata['source_sha256'] or
                    current['library_sha256'] != executor.build_metadata['library_sha256']):
                raise NativeBackendUnavailable('Cached native library changed; restart this process')
            _EXECUTORS.move_to_end(key)
    except NativeBackendUnavailable as exc:
        if diagnostics is not None:
            diagnostics.update(native_available=False, native_checked_points=0,
                               native_unavailable_reason=str(exc))
        return None
    if diagnostics is not None:
        diagnostics['native_available'] = True
    return executor.batch_density(alpha, invariants, m2, diagnostics=diagnostics)


def native_batch_density(compiled, alpha, invariants=None, m2=1, *, diagnostics=None):
    """Convenience native batch with direct standard-library Decimal fallback."""
    result=try_native_batch_density(compiled,alpha,invariants,m2,diagnostics=diagnostics)
    if result is None:
        # This entry point may itself be called by CompiledBPH.batch_density's
        # flagged subset. Calling that public batch again would recurse when the
        # optional library is absent. The standard-library Decimal path is the
        # direct, dependency-free fallback for this already-selected subset.
        result=[];highest=0
        for row in alpha:
            value,digits=compiled._adaptive_decimal(row,invariants,m2)
            result.append(float(value));highest=max(highest,digits)
        if diagnostics is not None:
            diagnostics.update(decimal_fallback_points=len(result),maximum_decimal_digits=highest)
        return result
    return result

"""Compare blastermaster output files with control files, numerically.

Text files are compared line by line: the numbers in each line must be within tolerance
and the text around them must be identical (ignoring whitespace, so that a number's width
or a line ending may change). The binary grids are compared as arrays of floats.
"""
import re
from dataclasses import dataclass
from typing import Optional

import numpy as np

NUMBER = re.compile(rb"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")
WHITESPACE = re.compile(rb"\s+")

# big-endian Fortran unformatted files: indices of the records holding float32 values
# (DelPhi phi maps: label, title, grid, label, scale & center)
FLOAT_RECORDS = {"vdw.vdw": (0, 1), "trim.electrostatics.phi": (2, 4)}


@dataclass
class Tolerance:
    rtol: float = 0.0
    atol: float = 0.0


@dataclass
class Result:
    file_name: str
    num_values: int = 0
    num_differing: int = 0  # values outside tolerance
    max_abs_diff: float = 0.0
    text_mismatch: Optional[str] = None  # first difference in the non-numeric content

    @property
    def ok(self):
        return self.num_differing == 0 and self.text_mismatch is None

    def __str__(self):
        s = f"{self.file_name}: {self.num_differing}/{self.num_values} values differ (max abs diff {self.max_abs_diff:.3g})"
        return s + (f"; {self.text_mismatch}" if self.text_mismatch else "")


def read_fortran_records(path):
    """Records of a big-endian Fortran unformatted sequential file."""
    with open(path, "rb") as f:
        data = f.read()
    records = []
    i = 0
    while i < len(data):
        n = int.from_bytes(data[i: i + 4], "big")
        records.append(data[i + 4: i + 4 + n])
        i += n + 8
    return records


def _split_numbers(data):
    """(text with the numbers taken out, per line; the numbers) of a text file's content."""
    skeletons, numbers = [], []
    for line in data.splitlines():
        skeletons.append(WHITESPACE.sub(b"", NUMBER.sub(b"#", line)))
        numbers.extend(float(n) for n in NUMBER.findall(line))
    return skeletons, np.array(numbers)


def _parts(path, file_name):
    """(non-numeric parts, numbers) of a file."""
    if file_name in FLOAT_RECORDS:
        records = read_fortran_records(path)
        float_indices = FLOAT_RECORDS[file_name]
        texts = [r for i, r in enumerate(records) if i not in float_indices]
        floats = [np.frombuffer(records[i], dtype=">f4") for i in float_indices]
        return texts, np.concatenate(floats).astype(np.float64)
    with open(path, "rb") as f:
        return _split_numbers(f.read())


def compare_files(expected_path, actual_path, file_name, tolerance=Tolerance()):
    result = Result(file_name)
    expected_texts, expected_values = _parts(expected_path, file_name)
    actual_texts, actual_values = _parts(actual_path, file_name)

    for n, (e, a) in enumerate(zip(expected_texts, actual_texts)):
        if e != a:
            result.text_mismatch = f"part {n + 1}: expected {e[:80]!r}, got {a[:80]!r}"
            break
    else:
        if len(expected_texts) != len(actual_texts):
            result.text_mismatch = f"expected {len(expected_texts)} lines/records, got {len(actual_texts)}"

    result.num_values = len(expected_values)
    if len(expected_values) != len(actual_values):
        result.num_differing = max(len(expected_values), len(actual_values))
        result.text_mismatch = result.text_mismatch or f"expected {len(expected_values)} values, got {len(actual_values)}"
        return result
    close = np.isclose(actual_values, expected_values, rtol=tolerance.rtol, atol=tolerance.atol, equal_nan=True)
    result.num_differing = int(np.count_nonzero(~close))
    if len(expected_values):
        result.max_abs_diff = float(np.nanmax(np.abs(actual_values - expected_values)))
    return result

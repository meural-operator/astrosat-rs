"""Independent pure-Python FITS reader for verification (zero shared code with astrosat-fits)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple
import numpy as np


_FITS_TYPE_MAP: Dict[str, Tuple[str, int]] = {
    "L": ("|S1", 1),
    "B": (">u1", 1),
    "I": (">i2", 2),
    "J": (">i4", 4),
    "K": (">i8", 8),
    "E": (">f4", 4),
    "D": (">f8", 8),
    "A": ("|S1", 1),
}


def _parse_card_value(raw: str) -> Any:
    raw = raw.strip()
    if raw.startswith("'"):
        # FITS string literal; handle escaped '' inside string
        out = []
        i = 1
        while i < len(raw):
            if raw[i] == "'":
                if i + 1 < len(raw) and raw[i + 1] == "'":
                    out.append("'")
                    i += 2
                    continue
                break
            out.append(raw[i])
            i += 1
        return "".join(out).rstrip()
    # Strip trailing comment
    val_part = raw.split("/", 1)[0].strip()
    if val_part == "T":
        return True
    if val_part == "F":
        return False
    try:
        return int(val_part)
    except ValueError:
        pass
    try:
        return float(val_part.replace("D", "E"))
    except ValueError:
        return val_part


def _parse_tform(tform: str) -> Tuple[int, str]:
    tform = tform.strip()
    i = 0
    while i < len(tform) and tform[i].isdigit():
        i += 1
    repeat = int(tform[:i]) if i > 0 else 1
    code = tform[i:i + 1]
    return repeat, code


@dataclass
class Hdu:
    index: int
    header: Dict[str, Any]
    cards: List[str]
    data_bytes: memoryview

    @property
    def extname(self) -> str:
        return str(self.header.get("EXTNAME", "PRIMARY")).strip()

    @property
    def extver(self) -> int:
        return int(self.header.get("EXTVER", 1))

    @property
    def xtension(self) -> str:
        return str(self.header.get("XTENSION", "IMAGE")).strip()

    def image_array(self) -> Optional[np.ndarray]:
        naxis = int(self.header.get("NAXIS", 0))
        if naxis == 0 or self.xtension == "BINTABLE":
            return None
        bitpix = int(self.header.get("BITPIX", 8))
        dtype_map = {
            8: ">u1",
            16: ">i2",
            32: ">i4",
            64: ">i8",
            -32: ">f4",
            -64: ">f8",
        }
        dims = [int(self.header[f"NAXIS{i + 1}"]) for i in range(naxis)]
        count = 1
        for d in dims:
            count *= d
        arr = np.frombuffer(self.data_bytes, dtype=dtype_map[bitpix], count=count)
        arr = arr.reshape(tuple(reversed(dims)))
        bscale = float(self.header.get("BSCALE", 1.0))
        bzero = float(self.header.get("BZERO", 0.0))
        if bscale != 1.0 or bzero != 0.0:
            arr = arr.astype(np.float64) * bscale + bzero
        return arr

    def bintable_columns(self) -> List[Tuple[str, int, str, int]]:
        """Returns list of (name, repeat, code, byte_offset)."""
        tfields = int(self.header.get("TFIELDS", 0))
        cols = []
        offset = 0
        for i in range(1, tfields + 1):
            name = str(self.header.get(f"TTYPE{i}", f"COL{i}")).strip()
            tform = str(self.header[f"TFORM{i}"]).strip()
            repeat, code = _parse_tform(tform)
            if code in ("P", "Q"):
                size = 8 if code == "P" else 16
            elif code == "X":
                size = (repeat + 7) // 8
            else:
                _, elem_size = _FITS_TYPE_MAP[code]
                size = repeat * elem_size
            cols.append((name, repeat, code, offset))
            offset += size
        return cols

    def read_column(self, col_name: str) -> np.ndarray:
        if self.xtension != "BINTABLE":
            raise ValueError(f"HDU {self.index} ({self.extname}) is not BINTABLE")
        nrows = int(self.header.get("NAXIS2", 0))
        row_len = int(self.header.get("NAXIS1", 0))
        target = col_name.strip().upper()
        for idx, (name, repeat, code, offset) in enumerate(self.bintable_columns(), start=1):
            if name.upper() == target:
                if code in ("P", "Q", "X"):
                    raise NotImplementedError(f"Column {col_name!r} has unsupported TFORM code {code!r} in fits_reader")
                if nrows == 0:
                    if code == "A":
                        return np.empty(0, dtype=f"U{repeat}")
                    dt, _ = _FITS_TYPE_MAP[code]
                    return np.empty((0, repeat) if repeat > 1 else 0, dtype=dt)
                raw = np.frombuffer(self.data_bytes[: nrows * row_len], dtype=np.uint8).reshape(nrows, row_len)
                if code == "A":
                    sl = raw[:, offset : offset + repeat].tobytes()
                    return np.array(
                        [sl[r * repeat : (r + 1) * repeat].decode("ascii", errors="replace").rstrip() for r in range(nrows)]
                    )
                if code == "L":
                    sl = raw[:, offset : offset + repeat]
                    arr = (sl == ord("T")).copy()
                    return arr[:, 0] if repeat == 1 else arr
                dt, elem_size = _FITS_TYPE_MAP[code]
                sl = np.ascontiguousarray(raw[:, offset : offset + repeat * elem_size])
                arr = np.frombuffer(sl.tobytes(), dtype=dt)
                if repeat > 1:
                    arr = arr.reshape(nrows, repeat)
                tscal = float(self.header.get(f"TSCAL{idx}", 1.0))
                tzero = float(self.header.get(f"TZERO{idx}", 0.0))
                if tscal != 1.0 or tzero != 0.0:
                    arr = arr.astype(np.float64) * tscal + tzero
                return arr
        raise KeyError(f"Column {col_name!r} not found in {self.extname}")


def read_fits(path: str) -> List[Hdu]:
    with open(path, "rb") as f:
        buf = memoryview(f.read())
    pos = 0
    hdus: List[Hdu] = []
    idx = 0
    while pos < len(buf):
        cards: List[str] = []
        hdr: Dict[str, Any] = {}
        while pos < len(buf):
            block = bytes(buf[pos : pos + 2880])
            pos += 2880
            end_found = False
            for i in range(0, 2880, 80):
                card = block[i : i + 80].decode("ascii", errors="replace")
                cards.append(card)
                if card.startswith("END     "):
                    end_found = True
                    break
                key = card[:8].strip()
                if len(card) >= 10 and card[8:10] == "= ":
                    hdr[key] = _parse_card_value(card[10:])
            if end_found:
                break
        if not hdr:
            break
        bitpix = abs(int(hdr.get("BITPIX", 8)))
        naxis = int(hdr.get("NAXIS", 0))
        pcount = int(hdr.get("PCOUNT", 0))
        gcount = int(hdr.get("GCOUNT", 1))
        if naxis == 0:
            data_bytes = 0
        else:
            dims = [int(hdr.get(f"NAXIS{i + 1}", 0)) for i in range(naxis)]
            prod = 1
            for d in dims:
                prod *= d
            data_bytes = (bitpix // 8) * gcount * (pcount + prod)
        pad = (2880 - (data_bytes % 2880)) % 2880
        data_view = buf[pos : pos + data_bytes]
        pos += data_bytes + pad
        hdus.append(Hdu(index=idx, header=hdr, cards=cards, data_bytes=data_view))
        idx += 1
    return hdus

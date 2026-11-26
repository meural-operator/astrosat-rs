"""Inspect official AstroSat CZTI & LAXPC Level-2 reference products."""
from __future__ import annotations

import os
import numpy as np
from fits_reader import read_fits


def inspect_file(path: str) -> None:
    print(f"\n=== {os.path.basename(path)} ===")
    hdus = read_fits(path)
    for h in hdus:
        if h.xtension == "BINTABLE":
            nrows = h.header.get("NAXIS2", 0)
            cols = [c[0] for c in h.bintable_columns()]
            extra = []
            for k in ("EXPOSURE", "TSTART", "TSTOP", "MJDREFI", "MJDREFF", "TIMESYS", "TIMEREF", "RA_PNT", "DEC_PNT", "GAINAPP"):
                if k in h.header:
                    extra.append(f"{k}={h.header[k]}")
            print(f"  HDU {h.index} [{h.extname}:{h.extver}] BINTABLE rows={nrows} cols={cols[:12]}{'...' if len(cols)>12 else ''} | {' '.join(extra)}")
            if h.extname == "GTI" and nrows > 0:
                starts = h.read_column("START")
                stops = h.read_column("STOP")
                print(f"    GTI exposure={np.sum(stops - starts):.4f} s, intervals={len(starts)}")
        else:
            arr = h.image_array()
            shape = arr.shape if arr is not None else ()
            s = float(np.sum(arr)) if arr is not None and arr.size > 0 else 0.0
            extra = [f"{k}={h.header[k]}" for k in ("EXPOSURE", "TSTART", "TSTOP", "RA_PNT", "DEC_PNT", "TWIST", "CRVAL1", "CRVAL2") if k in h.header]
            print(f"  HDU {h.index} [{h.extname}:{h.extver}] IMAGE shape={shape} sum={s:.4g} | {' '.join(extra)}")


def main() -> None:
    czti_dir = "/mnt/d/ISRO/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0"
    for fn in [
        "AS1T01_112T01_9000000406czt_level2_quad_clean.evt",
        "AS1T01_112T01_9000000406cztM0_level2_bc.evt",
        "AS1T01_112T01_9000000406czt_level2.gti",
        "AS1T01_112T01_9000000406cztM0_level2.aspect_Q0",
        "AS1T01_112T01_9000000406cztM0_level2_quad_clean.dph",
        "AS1T01_112T01_9000000406cztM0_level2_quad_clean.dpi",
        "AS1T01_112T01_9000000406cztM0_level2_quad_clean.img_Q0",
    ]:
        p = os.path.join(czti_dir, fn)
        if os.path.exists(p):
            inspect_file(p)

    lxp_dir = "/mnt/d/ISRO/data/sample/laxpc_l2/20160107_G02_010T01_9000000252_level2/laxpc/01509"
    for unit in ("lxp1", "lxp2", "lxp3"):
        for fn in [
            f"AS1G02_010T01_9000000252{unit}EA_level2.fits",
            f"AS1G02_010T01_9000000252{unit}EA_level2.lc",
            f"AS1G02_010T01_9000000252{unit}EA_level2.pha",
            f"AS1G02_010T01_9000000252{unit}EA_level2.gti",
        ]:
            p = os.path.join(lxp_dir, unit, "modeEA", fn)
            if os.path.exists(p):
                inspect_file(p)


if __name__ == "__main__":
    main()

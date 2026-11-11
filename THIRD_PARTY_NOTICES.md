# Third-Party Notices

zhyd is MIT-licensed (see `LICENSE`). This file lists the external work it
builds on, what of it (if anything) is included, and under which terms.

## reduce

- **Used in**: hydrogen placement, scoring heuristics and bond topology tables
  (`src/place/`, `src/optimize/`)
- **Source**: https://github.com/rlabduke/reduce
- **What is included**: no source code. zhyd is an independent
  reimplementation of the published method; the algorithm and parameters are
  owed to J. Michael Word and contributors (Word et al., J. Mol. Biol. 285,
  1735–1747, 1999).

## PROPKA 3 model

- **Used in**: the `empirical` pKa engine (`src/pka_empirical/`)
- **Source**: https://github.com/jensengroup/propka (PROPKA 3.5.1)
- **What is included**: no source code. The engine was written from a
  specification of the method (`docs/specs/zhyd-pka-empirical-engine.md`),
  the PROPKA 3 papers (Olsson et al., J. Chem. Theory Comput. 7, 525–537,
  2011; Søndergaard et al., J. Chem. Theory Comput. 7, 2284–2295, 2011) and
  black-box output of PROPKA, under the clean-room protocol in that
  specification. Numerical parameter values of the PROPKA 3.5.1 model are
  reproduced in parity mode and cited there.

## PROPKA 3.5.1 via zpropka (optional, LGPL-2.1)

- **Used in**: `--engine propka`, only when built with `zig build -Dpropka=true`
- **Source**: https://github.com/N283T/zpropka (a Zig port of PROPKA 3.5.1)
- **License**: LGPL-2.1-or-later, as PROPKA
- **What is included**: nothing in the default build. The dependency is lazy
  and fetched only with `-Dpropka=true`; binaries built that way are subject
  to the LGPL-2.1 for the zpropka part. Prebuilt or packaged zhyd
  distributions must not bundle it.

## zsasa

- **Used in**: patterns and code adapted where marked "Adapted from zsasa"
- **Source**: https://github.com/N283T/zsasa
- **License**: MIT, same author (Copyright (c) 2025 N283T)

## Chemical Component Dictionary (CCD)

- **Used in**: optional `-d components.cif` input for non-standard residues
- **Source**: wwPDB, https://www.wwpdb.org/data/ccd
- **What is included**: nothing; users supply the file. wwPDB data are
  available under CC0 1.0.

## Experimental pKa data used for parameter fitting

- **Used in**: the recalibrated zhyd parameter tables
  (`src/pka_empirical/params.zig`; fit records in `scripts/pka/params/`)
- **Sources**: PKAD-R (Chen et al., J. Comput. Biophys. Chem. 24, 1189–1203,
  2025; https://compbio.clemson.edu/pkad-r/) and PKAD-3_V2 with the
  time-stamped set from KaML-ESM (Shen et al.; CC BY-NC 4.0,
  https://github.com/JanaShenLab/KaML-ESM)
- **What is included**: no data. Only a few fitted parameter values and
  summary statistics are distributed; the scripts in `scripts/pka/` rebuild
  the data sets from the original downloads.

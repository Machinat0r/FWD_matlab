# MMS KH FOTE / FOTE-V

`kh_fote_event.m` applies the same four-point affine expansion to the magnetic
field and to a plasma velocity field.  The default velocity is the ion bulk
velocity `Vi`; use `VelocityField='Ve'` to switch to the electron bulk velocity.

The default is **no smoothing** (`SmoothSamples=1`).  Values greater than one
apply one moving-median pass only.  Time interpolation is limited to contiguous
data segments and never extrapolates across burst gaps.

The preview command is:

```matlab
result = run_kh005_preview;
```

Magnetic FOTE quality uses both `eta=|div B|/|curl B|` and the eigenvalue-trace
residual `xi`.  FOTE-V does not force `div(V)=0`; it uses the continuity proxy
`alpha=|div(nV)|/|curl(nV)|`.  Null distances are normalized by the mean of the
six tetrahedron edge lengths `L`, and the preview defines `d_min <= 2L` as near.

The generic call is:

```matlab
result = kh_fote_event(eventID,startUtc,endUtc, ...
    'VelocityField','Vi','SmoothSamples',1);
```

`kh_fote_event.py` is a numerically equivalent direct-CDF runner.  It is used
on this workstation because MATLAB R2026a currently crashes inside the IRFU
`EpochTT`/CDF bridge.  The Python runner keeps the same four-point expansion,
type criteria, gap policy, quality limits and default no-smoothing policy:

```powershell
python kh_fote_event.py KH005 2015-10-01T18:01:24Z 2015-10-01T18:09:00Z
```

For a physical-time average, use `--smooth-seconds`.  Magnetic field and ion
velocity are then averaged separately on their native timestamps in an exact
centered window before both products are aligned to the common ion timeline.
These runs add two panels at the top containing the unsmoothed magnetic field
and ion velocity, followed by the time-averaged fields and the FOTE/FOTE-V
results.  Finite interpolation is restricted to the actual source support, so
the moving average and plotting never extend across burst gaps.

The Python runner also ports the Poincare-index calculation from IRFU-MATLAB
`c_4_poincare_index.m` and `irf/+irf/solidangle.m`.  It preserves the IRFU
oriented-face order 123, 142, 134, and 243.  The same topological-degree
calculation is applied to the four-spacecraft magnetic and ion-velocity
vectors and shown as `PI_B` and `PI_V` in one panel.  Values near +1 or -1
identify an odd number of enclosed 3-D zeros; zero means zero or an even
number.  Both indices are exported to the CSV and MAT files.
For KH005, a 10 s window corresponds to about 1280 FGM samples at 128 Hz and
67 DIS samples at 6.667 Hz:

```powershell
python kh_fote_event.py KH005 2015-10-01T18:01:24Z 2015-10-01T18:09:00Z --smooth-samples 1 --smooth-seconds 10
```

The corresponding 5 s test uses about 640 FGM samples and 33 DIS samples:

```powershell
python kh_fote_event.py KH005 2015-10-01T18:01:24Z 2015-10-01T18:09:00Z --smooth-samples 1 --smooth-seconds 5
```

The current quality-only batch layout uses the FOTE/FOTE-V 40% error criteria,
does not reject or accept a result based on null distance, replaces the two
distance panels with magnetic FOTE (eta and xi) and FOTE-V (alpha) error
panels, keeps only the event ID as the title, and writes a multi-page PDF.
Use one sample for an unsmoothed run:

```powershell
python kh_fote_batch.py --smooth-samples 1
```

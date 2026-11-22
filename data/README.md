# Data guide

| File | Content |
|---|---|
| [classes.json](classes.json) · [readable table](all_classes.md) | Primitive trace 3–60 Sol conjugacy classes with monodromies, reversal witnesses, eta, signed minima and oriented orbit representatives. |
| [thresholds.json](thresholds.json) | Four full-line Sol certificates with Hermite coefficients, moments, comparison bounds, origin and dimension of the first signed shells. |
| [shells.json](shells.json) | Fifteen spherical spin lifts with rational coefficients of the multiplicity difference on each congruence class, ordinary eta and the two contact multiplicities. |
| [lens.csv](lens.csv) | All 434 lifted round lens actions through group order 32 in the rotation-weight convention. |
| [flat.json](flat.json) · [readable table](flat_thresholds.md) | Six orientable Bieberbach types, 28 lifts and critical aspect ratios. |
| [nil.json](nil.json) · [readable table](nil_thresholds.md) | 24 spin structures on the four Heisenberg quotients evaluated at three aspect ratios. Includes projected and full thresholds. |
| [nil_measurements.json](nil_measurements.json) | LC- and left/right-detector occupations for specified masses, intervals, modes, charges and Holst coefficients. |
| [device_certificate.json](device_certificate.json) | Four ball enclosures of positive/negative-velocity comparisons of the representative shell and ratios. Positive velocity means negative NY charge. |

Every file here, including the two readable tables, is written by `python -B script/paper08.py recompute`; `verify` checks the distributed copies, including that the readable tables agree with their JSON sources. Floating-point diagnostics in nil_measurements.json may differ in the last digits between LAPACK builds. [Reproduction](../script/README.md)

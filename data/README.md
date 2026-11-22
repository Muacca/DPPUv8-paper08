# Data guide

| File | Content |
|---|---|
| [classes.json](classes.json) · [readable table](all_classes.md) | Primitive trace 3–60 Sol conjugacy classes with monodromies, reversal witnesses, eta, signed minima and oriented orbit representatives. |
| [thresholds.json](thresholds.json) | Four full-line Sol certificates with Hermite coefficients, moments, comparison bounds, origin and dimension of the first signed shells. |
| [shells.json](shells.json) | Fifteen spherical spin lifts with rational coefficients of the multiplicity difference on each congruence class, ordinary eta and the two contact multiplicities. |
| [lens.csv](lens.csv) | All 434 lifted round lens actions through group order 32 in the rotation-weight convention. |
| [flat.json](flat.json) · [readable table](flat_thresholds.md) | Six orientable Bieberbach types, 28 lifts and critical aspect ratios. The stored JSON additionally retains sample expansion-window comparisons from the source calculation; only its spins table is used or regenerated here. |
| [nil.json](nil.json) · [readable table](nil_thresholds.md) | 24 spin structures on the four Heisenberg quotients evaluated at three aspect ratios. Includes projected and full thresholds. |
| [nil_measurements.json](nil_measurements.json) | LC- and left/right-detector occupations for specified masses, intervals, modes, charges and Holst coefficients. |
| [device_certificate.json](device_certificate.json) | Four ball enclosures of positive/negative-velocity comparisons of the representative shell and ratios. Positive velocity means negative NY charge. |

Source calculations are the project analyses on the Sol arithmetic census (2026-11-22), unshifted spherical quotient character tables (2026-11-21), first lens-shell zeros (2026-11-21), flat and Nil quotient surveys and charged Nil transport (2026-11-22). All data needed to run the public supplement are local. [Reproduction](../script/README.md)

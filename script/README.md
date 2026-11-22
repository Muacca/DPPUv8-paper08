# Reproducing Paper 08

From the repository root, use Python 3.11 or later and install the scientific packages in [requirements.txt](requirements.txt). Tested with Python 3.12.10, NumPy 2.3.5, SciPy 1.16.3, python-flint 0.8.0 and Matplotlib 3.10.7.

~~~sh
python -m pip install -r script/requirements.txt
python -B script/paper08.py verify
python -B script/paper08.py recompute --output results_new
python -B script/paper08.py render --data results_new --output figures_new
~~~

The recompute and render commands each take a new output directory. The render command reads only the selected data. Its PDF is the two-panel figure in the manuscript. Recompute enumerates the Sol classes, evaluates the Hermite residual certificates, computes the spherical shell tables and lens-energy zeros, enumerates the flat and Nil reference spectra and integrates the massive-shell experiment. The interval certificate uses python-flint's Arb arithmetic alongside the numerical SciPy propagator. Numerical representative vectors may differ by a sign or rounding between LAPACK builds. The physical conclusions and interval isolation, not byte identity of floating-point intermediate files, are checked.

The verify command reads the supplied data, checks all 488 stored integer representatives and eta conventions, the number of spin lifts, signed threshold examples, and inclusion of the detector values in their interval certificates. The English article provides the operator argument needed for arbitrary Sol trace. See the [data scope](../data/SCOPE.md).

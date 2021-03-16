# CXRO_calc
This repo uses physical properties from the [online CXRO atomic scattering factor database](https://henke.lbl.gov/optical_constants/) to calculate properties that were useful for my doctoral research. The structure of this repo is as follows:

- Local copies of CXRO database are located in the [CrossSections](CrossSections), [index](index), and [Transmissions](Transmissions) directories.
- Convenience functions are located in the [funcs](funcs) directory.
- Plots are generated using the scripts in the [test](test) directory, with results saved in the [figures](figures) directory.

The scripts are described below:
- [CXRO.py](test/CXRO.py): transmission of extreme ultraviolet (XUV) light through gas media, metallic filters and samples
- [critical_fraction.py](test/critical_fraction.py): critical ionization and phase matching calculations for high harmonic generation (HHG)
- [Fresnel.py](test/Fresnel.py): reflectivity of mirrors in the XUV range as a function of material, incidence angle, and mirror smoothness
- [EM_geometry.py](test/EM_geometry.py): to-scale drawings of our XUV mirror and its deviation from an ideal ellipsoid
- [real_imag_index_plotting.py](test/real_imag_index_plotting.py): plots the relative error by ignoring edge-effects when computing XUV transmission in a thin film sample

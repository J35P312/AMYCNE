# AMYCNE
========
AMYCNE is a copy number estimation toolkit, designed for WGS data. It contains three modules. genotyping, which is used to estimate the copy number in one or more target region. annotation, which is used to estimate the copy number across structural variants stored in a vcf file, and annotates each variant with its copy number. The third module is the call module, this module is used to detect structural variants. The call module is still in early development.
# Installation
========
AMYCNE requires scipy and numpy. AMYCNE has been tested on python 2.7.11, byt might run on older versions of python as well.
To improve the performance of AMYCNE, the code of AMYCNE may be compiled using cython:

python setup.py build_ext --inplace

# Input
========
For each mode, AMYCNE requires at least two input file. One input file describing the GC in bins across the entire genome, and one file describing the coverage in bins of the same size as those in the GC content file. The GC content file may be created using the Generate_GC_tab.py, but any tab separated GC file should work.
The coverage file could be generated using TIDDIT, which sambaba might work as well.


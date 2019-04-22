
# Auxiliary scripts for Quanformer
Last updated: Apr 22 2019

| Script                | Description                                                                            |
| ----------------------|----------------------------------------------------------------------------------------|
| `cat_mols.py`         | concatenates molecules from various files into single output file, modified from [OpenEye](https://docs.eyesopen.com/toolkits/python/_downloads/catmols.py) |
| `mols2pdf.py`         | draws molecules for PDF output, modified from [OpenEye](https://docs.eyesopen.com/toolkits/python/_downloads/mols2pdf.py) |
| `run_local.sh`        | bash script to iterate over Psi4 calculations on a local machine                       |
| `viewer.ipynb`        | visualize molecules in iPython notebook                                                |


## TO BE TRANSFERRED
| Script                | Description                                                                            |
| ----------------------|----------------------------------------------------------------------------------------|
| `align_two_mols.py`   | read in a molecule from two distinct files and match their atom indices                |
| `cleanfromvmd.py`     | clean molecules that were filtered ad hoc through VMD                                  |
| `convertExtension.py` | read in mol from one file type and write out to another file type                      |
| `findIntraHB.py`      | identify molecules that may have internal hydrogen bonding                             |
| `jobcount.sh`         | count total number of conformer calculations as well as number of unfinished jobs      |
| `loadFromXYZ.py`      | copy coordinates from XYZ to MOL2 file                                                 |
| `selectConfs.tcl`     | script for VMD to further filter molecule set, e.g., by some internal distance         |
| `write_first_confs.py`| write out first conformers of all molecules                                            |
| `writeOneMol.py`      | write out single mol and all its conformers OR single conformer of specified mol       |
| `xyzByStep.sh`        | simple Bash processing of Psi4 output file to see geometries throughout optimization   |

Other helpful scripts:
 * [This list](http://vergil.chemistry.gatech.edu/resources/utilities.html) from the Sherrill group
 * [OpenEye script](https://docs.eyesopen.com/toolkits/python/oechemtk/oechem_examples/oechem_example_molextract.html) to extract molecules by title


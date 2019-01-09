
# Data sources for testing
Notes for VTL.  
Last updated: Dec 19 2018

| filename                       | test script          | source                                                                                       |
|--------------------------------|----------------------|----------------------------------------------------------------------------------------------|
| `carbon-222.sdf`,`t1`,`s1`     | `get_psi_results.py` | `/beegfs/DATA/mobley/limvt/openforcefield/pipeline/work_hydrocarbons/HESS`                   |
| `cooh`                         | `getTurbResults.py`  | `/beegfs/DATA/mobley/limvt/openforcefield/pipeline/work_coohDirection/02_torsion/03_Turbomole/2_solv-HF` |
| `freeze.sdf`                   | `confs_to_psi.py`    | `/data11/home/jmaat/off_nitrogens/sdf_min/sdf_min_mol2/pyrnit_2_constituent_11_improper.sdf` |
| `GBI`                          | `get_psi_results.py` | `/beegfs/DATA/mobley/limvt/openforcefield/pipeline/03_examples/set1/GBI`                     |
| `gbi-200.sdf`,`gbi_single.sdf` | `get_psi_results.py` | `/beegfs/DATA/mobley/limvt/openforcefield/pipeline/03_examples/set1/examples2-200.sdf`       |
| `gbi_single.sdf`               | `initialize_confs.py`| see above (same file)                                                                        |
| `gbi.sdf`                      | `filter_confs.py`    | `/beegfs/DATA/mobley/limvt/openforcefield/pipeline/03_examples/set1/examples2.sdf`           |
| `methane.smi`                  | `initialize_confs.py`| self-generated                                                                               |
| `methane_c2p.sdf`              | `confs_to_psi.py`    | `test_initialize_confs()`                                                                           |
| `methane_title-1.0.sdf`        | `get_psi_results.py` | edited copy of `methane_c2p.sdf`                                                             |
| `output_hess.dat`              | `get_psi_results.py` | `/beegfs/DATA/mobley/limvt/openforcefield/hessian/sandbox_benzene/benzene/output.dat`        |
| `output_opt.dat`, `timer.dat`  | `get_psi_results.py` | `/beegfs/DATA/mobley/limvt/openforcefield/pipeline/03_examples/set1/GBI/1/`                  |
| `output_spe.dat`               | `get_psi_results.py` | `/beegfs/DATA/mobley/limvt/openforcefield/pipeline/set1_01_main/SPE2/AlkEthOH_c1178/1/output.dat` |
| `steric_clash.smi`             | `initialize_confs.py`| `/DFS-L/old_beegfs_data/mobley/limvt/openforcefield/pipeline/set1_01_main/set1_01_main.smi`  |


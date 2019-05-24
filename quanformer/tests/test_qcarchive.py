import pytest

import qcfractal.interface as ptl
from qcfractal.interface.models import GridOptimizationInput
from qcfractal.testing import fractal_compute_server, using_geometric, using_rdkit

@using_geometric
@using_rdkit
def test_service_gridoptimization_single_noopt(fractal_compute_server):

    client = ptl.FractalClient(fractal_compute_server)

    # Add a HOOH
    hooh = ptl.data.get_molecule("hooh.json")
    initial_distance = hooh.measure([1, 2])

    # Options
    service = GridOptimizationInput(**{
        "keywords": {
            "preoptimization": False,
            "scans": [{
                "type": "distance",
                "indices": [1, 2],
                "steps": [-0.1, 0.0],
                "step_type": "relative"
            }]
        },
        "optimization_spec": {
            "program": "geometric",
            "keywords": {
                "coordsys": "tric",
            }
        },
        "qc_spec": {
            "driver": "gradient",
            "method": "UFF",
            "basis": "",
            "keywords": None,
            "program": "rdkit",
        },
        "initial_molecule": hooh,
    }) # yapf: disable

    ret = client.add_service([service])
    fractal_compute_server.await_services()
    assert len(fractal_compute_server.list_current_tasks()) == 0

    result = client.query_procedures(id=ret.ids)[0]

    assert result.status == "COMPLETE"
    assert result.starting_grid == (1, )
    assert pytest.approx(result.get_final_energies((0, )), abs=1.e-4) == 0.00032145876568280524

    assert result.starting_molecule == result.initial_molecule

    # Check initial vs startin molecule
    assert result.initial_molecule == result.starting_molecule

    mol = client.query_molecules(id=result.starting_molecule)[0]
    assert pytest.approx(mol.measure([1, 2])) == initial_distance

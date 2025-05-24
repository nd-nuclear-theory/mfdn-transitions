""" runmfdn07.py

    See runmfdn.txt for description.

    Mark A. Caprio
    University of Notre Dame

    - 06/08/17 (pjf): Created, copied from runmfd01; switch to MFDn v15.
    - 07/31/17 (pjf): Set MFDn driver module in task dictionary.
    - 08/11/17 (pjf): Update for split single-particle and many-body truncation modes.
    - 09/24/17 (pjf): Save wavefunctions.
    - 12/19/17 (pjf): Update for mfdn->ncci rename.
    - 09/07/19 (pjf): Remove Nv from truncation_parameters.
    - 09/11/19 (pjf): Fix task-data save.
"""

import mcscript
import ncci
import ncci.mfdn_v15
import ncci.postprocessing

# initialize mcscript
mcscript.init()

##################################################################
# build task list
##################################################################

ncci.environ.interaction_run_list = [
    "example-data",
    "run0164-JISP16-ob-9",
    "run0164-JISP16-ob-13",
    "run0164-JISP16-tb-10",
    "run0164-JISP16-tb-20",
    "run0306-N2LOopt500",  # up to tb-20
    "runmac0471-Daejeon16-tb-20",
    "runvc0083-Daejeon16-ob-13"
]

task = {
    # nuclide parameters
    "nuclide": (2, 1),

    # Hamiltonian parameters
    "interaction": "Daejeon16",
    "use_coulomb": True,
    "a_cm": 0.,
    "hw_cm": None,

    # input TBME parameters
    "truncation_int": ("ob", 13),
    "hw_int": 20.,
    "truncation_coul": ("ob", 13),
    "hw_coul": 20.,

    # basis parameters
    "basis_mode": ncci.modes.BasisMode.kDirect,
    "hw": 20.,

    # transformation parameters
    "xform_truncation_int": None,
    "xform_truncation_coul": None,
    "hw_coul_rescaled": None,
    "target_truncation": None,

    # traditional oscillator many-body truncation
    "sp_truncation_mode": ncci.modes.SingleParticleTruncationMode.kNmax,
    "mb_truncation_mode": ncci.modes.ManyBodyTruncationMode.kFCI,
    "truncation_parameters": {
        "Nmax_orb": 4,
        "M": 0.5,
        "Nmax": 4,
        "Nstep": 2,
        },

    # diagonalization parameters
    "diagonalization": True,
    "eigenvectors": 5,
    "initial_vector": -2,
    "max_iterations": 200,
    "tolerance": 1e-6,
    "partition_filename": None,

    # obdme parameters
    ## "hw_for_trans": 20,
    "save_obdme": False,

    # two-body observables
    ## "tb_observable_sets": ["H-components","am-sqr"],
    "tb_observable_sets": ["H-components", "am-sqr", "isospin"],
    "tb_observables": [("Ntotal",  (0,0,0), ncci.operators.tb.Ntotal(A=3, hw=20.))],

    # version parameters
    "h2_format": 15099,
    "mfdn_executable": "v15-beta02/xmfdn-h2-lan",
    "mfdn_driver": ncci.mfdn_v15,

}

################################################################
# run control
################################################################

# add task descriptor metadata field (needed for filenames)
task["metadata"] = {
    "descriptor": ncci.descriptors.task_descriptor_7(task)
    }

#ncci.radial.set_up_interaction_orbitals(task)
#ncci.radial.set_up_orbitals(task)
#ncci.radial.set_up_xforms_analytic(task)
#ncci.radial.set_up_obme_analytic(task)
#ncci.tbme.generate_tbme(task)
#ncci.mfdn_v15.run_mfdn(task)
#ncci.mfdn_v15.save_mfdn_task_data(task)
#ncci.postprocessing.evaluate_ob_observables(task)
ncci.handlers.task_handler_oscillator_pre(task)
ncci.handlers.task_handler_oscillator_mfdn(task)
ncci.handlers.task_handler_post_run(task, cleanup=False)

##################################################################
# task control
##################################################################

## mcscript.task.init(
##     tasks,
##     task_descriptor=ncci.descriptors.task_descriptor_7,
##     task_pool=task_pool,
##     phase_handler_list=[ncci.handlers.task_handler_oscillator]
##     )

################################################################
# termination
################################################################

mcscript.termination()

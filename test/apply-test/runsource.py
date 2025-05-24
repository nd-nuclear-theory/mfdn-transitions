"""runsource.py

    qsubm --here source --pool=ALL --phase=0 --serialthreads=8 --threads=8
    qsubm --here source --pool=ALL --phase=1 --serialthreads=8 --threads=8

    Patrick J. Fasano
    University of Notre Dame

    - Based on runmfdn10.  Committed 05/10/21 (pjf). 
    - 05/21/25 (mac):
      - Update handler api.
      - Retreat from Nmax04/06 to Nmax02/04, to use example tb-6 interaction files.

    Would it be more transparent to split source and target generation into two
    different tasks with appropriate descriptors.

    Need to force same sp orbitals?

"""

import mcscript
import mcscript.ncci as ncci

# initialize mcscript
mcscript.init()

##################################################################
# build task list
##################################################################

ncci.environ.interaction_run_list = [
    "jisp16-tb-6",
    "coulomb-tb-6",
]

task = {
    # nuclide parameters
    "nuclide": (3, 3),

    # Hamiltonian parameters
    "interaction": "JISP16",
    "use_coulomb": True,
    "a_cm": 20.,
    "hw_cm": None,

    # input TBME parameters
    "truncation_int": ("tb", 6),
    "hw_int": 20.,
    "truncation_coul": ("tb", 6),
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
    "mb_truncation_mode": ncci.modes.ManyBodyTruncationMode.kNmax,
    "truncation_parameters": {
        "Nmax": 4,
        "Nstep": 2,
        "M": 0,
        },

    # diagonalization parameters
    "diagonalization": True,
    "eigenvectors": 2,
    "initial_vector": -2,
    "max_iterations": 200,
    "tolerance": 1e-6,
    "partition_filename": None,
    ## "save_tbme": True,
    
    # obdme parameters
    ## "hw_for_trans": 20,
    "obdme_multipolarity": 2,
    # "obdme_reference_state_list": [(0, 0, 1)],
    "calculate_obdme": False,
    "save_obdme": False,

    # two-body observables
    ## "tb_observable_sets": ["H-components","am-sqr"],
    "tb_observable_sets": ["H-components", "am-sqr"],
    "tb_observables": [
        ("J", (1,0,0), {"U[j]": 1.0}),
        ("identity", (0,0,0), {"identity": 1.0}),
    ],

    # wavefunction storage
    "save_wavefunctions": True,

    # version parameters
    "h2_format": 15099,
    "mfdn_executable": "xmfdn-h2-lan",
    "mfdn_driver": ncci.mfdn_v15
}

def handler0(task):
    """ Generate source wave functions (Nmax=4), with postfix 0."""
    ## ncci.handlers.task_handler_oscillator(task, postfix="0", cleanup=False)
    ncci.handlers.task_handler_mfdn_pre(task, postfix="0")
    ncci.handlers.task_handler_mfdn_run(task, postfix="0")
    ncci.handlers.task_handler_mfdn_post_no_cleanup(task, postfix="0")  # nocleanup is to save TBME file for operator

def handler1(task):
    """ Generate target wave function indexing (Nmax=6), with postfix 1."""
    task["truncation_parameters"]["Nmax"] += 2  # 05/21/25 (mac): Doesn't this leave us with a mismatch of sp orbitals with respect to the Nmax=4 indexing?
    task["h2_format"] = 15200
    ncci.handlers.task_handler_mfdn_pre(task, postfix="1")
    ncci.handlers.task_handler_mfdn_dimension(task, postfix="1")

mcscript.task.init(
    [task],
    task_descriptor=lambda task: "task",
    phase_handler_list=[handler0,handler1]
    )

################################################################
# termination
################################################################

mcscript.termination()

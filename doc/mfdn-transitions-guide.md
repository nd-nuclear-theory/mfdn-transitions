# Guide to `mfdn-transitions`

+ 07/16/21 (mac): Created, with notes on HPC run parameters.

+ 06/27/25 (mac): Add notes on executables and namelist input files.

----------------------------------------------------------------

## Executables

There are two basic executables:

   - `xtransitions`: Evaluates one-body densities and matrix elements of
     two-body operators, between two wave functions.
   
   - `xapply`: Applies a two-body operator to a wave function, to obtain a new
     wave function.
     
These basic executables are defined for Tz-conserving operators, that is, the
initial and final wave functions have the same Tz.  Both executables then have a
`-deltaTz` variant, which is defined to work with initial and final wave
functions differing by 1-2 units in Tz.
   

## Namelist input parameters

Run parameters are provided to the code through a Fortran namelist input file,
either named `transitions.input` (for `xtransitions`) or `apply.input` (for
`xapply`).  Example input files are provided in `doc/examples`.  Here we run
through an example of each type of input file, to highlight the meanings of the
parameters.

The full set of available parameters for `transitions.input` is defined as
follows (in `module_Transitions_input.f90`):

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Basis parameters
  real(4) :: fmass=938.92, hbomeg=0.0
  !
  ! Number of transition operators
  integer :: numTBtrans=0
  !
  ! Transition OBDMEs
  logical :: obdme=.true.
  integer :: max2K=4
  !
  ! <bra| and |ket>  state
  integer :: TwoJ_bra
  integer :: n_bra=1
  integer, dimension(MaxNumKets) :: TwoJ_ket=-1
  integer, dimension(MaxNumKets) :: n_ket=0
  !
  character(LEN=MAX_PATHLENGTH) :: infofilename_bra='mfdn_smwf.info'
  character(LEN=MAX_PATHLENGTH) :: infofilename_ket='mfdn_smwf.info'
  character(LEN=MAX_PATHLENGTH) :: basisfilename_bra='mfdn_MBgroups'
  character(LEN=MAX_PATHLENGTH) :: basisfilename_ket='mfdn_MBgroups'
  character(LEN=MAX_PATHLENGTH) :: smwffilename_bra='mfdn_smwf'
  character(LEN=MAX_PATHLENGTH) :: smwffilename_ket='mfdn_smwf'
  !
  character(LEN=MAX_PATHLENGTH), dimension(MaxNumTBMEops) :: TBMEoperators
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us take an example `transitions.input` file:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  &transition_data
  
  ! bra
  infofilename_bra = '../data/Z3-N3-Daejeon16-coul1-hw15.000-a_cm50-Nmax02-Mj1.0-lan600-tol1.0e-06/mfdn_smwf.info',
  basisfilename_bra = '../data/Z3-N3-Daejeon16-coul1-hw15.000-a_cm50-Nmax02-Mj1.0-lan600-tol1.0e-06/mfdn_MBgroups',
  smwffilename_bra = '../data/Z3-N3-Daejeon16-coul1-hw15.000-a_cm50-Nmax02-Mj1.0-lan600-tol1.0e-06/mfdn_smwf',
  TwoJ_bra = 2,
  n_bra = 1,
  
  ! ket
  infofilename_ket = '../data/Z3-N3-Daejeon16-coul1-hw15.000-a_cm50-Nmax02-Mj1.0-lan600-tol1.0e-06/mfdn_smwf.info',
  basisfilename_ket = '../data/Z3-N3-Daejeon16-coul1-hw15.000-a_cm50-Nmax02-Mj1.0-lan600-tol1.0e-06/mfdn_MBgroups',
  smwffilename_ket = '../data/Z3-N3-Daejeon16-coul1-hw15.000-a_cm50-Nmax02-Mj1.0-lan600-tol1.0e-06/mfdn_smwf',
  TwoJ_ket(1) = 2,
  n_ket(1) = 1,
  
  ! ob densities
  obdme = .true.,
  max2K = 4,
  hbomeg = 15,
  
  ! tb operators
  numTBtrans = 2,
  TBMEoperators(1) = '../data/identity_h2v15099',
  TBMEoperators(2) = '../data/Tintr_h2v15099',
  /
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `bra`: The parameters in this section specify the bra state to be used in
  computing the densities and matrix elements.
  
  + `infofilename_bra`, `basisfilename_bra`, `smwffilename_bra`: The paths to
    the bra wave function "info", "basis", and "smwf" (that is, "shell model
    wave function") files, respectively.  These files will most commonly all be
    found in the same directory, and with the standard filenames given to them
    when MFDn writes out its wave functions.  These are, namely,
    `mfdn_smwf.info`, `mfdn_MBgroups<nnnn>`, and `mfdn_smwf<nnnn>`,
    respectively.  Here `<nnnn>` represents the wave function segment index, as
    defined in, e.g., Fig. 1 of [aktulga2014:mfdn-scalability].  Only the base
    filename, before the segment index, should be specified in
    `basisfilename_bra` or `smwffilename_bra`.

  + `TwoJ_bra`, `n_bra`: The (J,n) quantum numbers for the bra state to select
    from the wave function files.  Here `TwoJ_bra` specifies 2*J, while `n` is
    the 1-based index.  Thus, e.g., in this example the first 1+ state ("1^+_1"
    in spectroscopic notation), is specified.
    
- `ket`: The parameters in this section similarly the ket state(s) to be used in
  computing the densities and matrix elements.
  
  + `infofilename_ket`, `basisfilename_ket`, `smwffilename_ket`: The paths to
    the ket wave function files (defined similarly to the bra parameters above).
  
  + `TwoJ_ket(i)`, `n_ket(i)`: The ket state(s) to be selected.  However, an
    important difference arises relative to the analogous parameters for the bra
    state above, due to the "multi-ket" computational capability of the
    postprocessor.  The postprocessor uses vectorization to simultaneously
    compute matrix elements from one bra state to multiple ket states, taken
    from the same wave function file, at no additional compute cost.  Presently,
    this functionality supports simultaneous calculations with up to 8 ket
    states, but the limit on the number of ket states is ultimately determined
    by the hardware vectorization capabilities.  Thus, the `TwoJ_ket(i)`
    `n_ket(i)` parameters are indexed, with the index running from 1 to 8.  In
    this example, only a single ket state, again, the first 1+ state, is
    specified.

- `ob densities`: The parameters in this section apply to the calculation of
  one-body densities.
  
  + `obdme`: Whether or not to compute densities.
  
  + `max2K`: Highest multipolarity k for the densities <|| (a+ a)_k ||> to be
    computed, given as 2*k.  In this example, densities through quadrupole (k=2)
    are requested.
    
  + `hbomeg`, `fmass`: Harmonic oscillator basis hw and mass parameters.  These
    are vestigal parameters from the oscillator-basis one-body observable
    calculation code in MFDn.  They do *not* affect the principal output
    (densities or two-body operator matrix elements) generated by the
    postprocessor in `transitions.res`, which makes no assumptions regarding the
    radial wave functions of the single-particle basis.  These parameters only
    affect some supplemental output provided in `transitions.out` (namely,
    multipole transition matrix elements generated assuming an oscillator
    basis).

- `tb operators`: The parameters in this section apply to the calculation of
  matrix elements of two-body operators.
  
  + `numTBtrans`: The number of two-body operators to consider.
  
  + `TBMEoperators(i)`: The paths to the two-body matrix element (TBME) files
  for these operators.  These should be h2 files, in either v15099 or v15200
  format.  The extension `.dat` or `.bin` should be omitted from the filename
  specified here in the namelist input.  The postprocessor will search for a
  file with either extension, and assume text or binary mode input, accordingly.

The available parameters for `apply.input` are similar to those for
`transitions.input`, but a few of those are unused, and a few new parameters are
defined as follows (in `MFDn_apply_operator.f90`):

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer, dimension(MaxNumKets) :: TwoJ_out=-1
  integer, dimension(MaxNumKets) :: n_out=0
  logical :: normalize=.false.
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us take an example `apply.input` file:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  &transition_data
  
  ! ket (source)
  infofilename_ket = '../data/Z3-N3-Daejeon16-coul1-hw15.000-a_cm50-Nmax02-Mj1.0-lan600-tol1.0e-06/mfdn_smwf.info',
  basisfilename_ket = '../data/Z3-N3-Daejeon16-coul1-hw15.000-a_cm50-Nmax02-Mj1.0-lan600-tol1.0e-06/mfdn_MBgroups',
  smwffilename_ket = '../data/Z3-N3-Daejeon16-coul1-hw15.000-a_cm50-Nmax02-Mj1.0-lan600-tol1.0e-06/mfdn_smwf',
  TwoJ_ket(1) = 2,
  n_ket(1) = 1,
  
  ! bra (target)
  infofilename_bra = 'target/mfdn_smwf.info',
  basisfilename_bra = 'target/mfdn_MBgroups',
  smwffilename_bra = 'target/mfdn_smwf',
  TwoJ_out(1) = 2,
  n_out(1) = 1,
  
  ! operator application
  TBMEoperators(1) = '../data/identity_h2v15099',
  normalize = .false.,
  
  /
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `ket (source)`: These parameters specify the source wave function(s), to which
  the two-body operator should be applied.  See description of `ket` parameters
  for `transitions.input` above.
  
- `bra (target)`: These parameters specify the target target wave functions.

  The interpretation of the files is a bit subtle, since varoius of these files
  serve as both inputs and output.  In preparation for running the
  postprocessor, one will typically create a target directory which contains a
  *template* "info" file `mfdn_smwf.info` and the appropriate "basis" files
  `mfdn_MBgroups<nnnn>`.
  
  The "info" file is first read as an input, to specify various truncation
  parameters defining the many-body basis for the output wave function.  But the
  list of state quantum numbers appearing at the end of this file is ignored.
  Then this file is *overwritten*, to contain a new state list at its end,
  reflecting the target wave functions actually written.  In the example
  considered here, the template file contains:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     15200    ! Version Number
     3   3   2    ! Z, N, 2Mj
  ...

     1      4.1000             711       1    ! Par, WTm, Dim, Npe

         8    ! n_states, followed by (i, 2J, nJ, T, -Eb, res)
         1       2       1    0.00        -25.5777        0.65E-05
         2       6       1    0.00        -24.1482        0.61E-05
         3       4       1    1.00        -19.9833        0.57E-05
         4       4       2    0.00        -18.9755        0.55E-05
         5       2       2    0.00        -17.0419        0.50E-05
         6       4       3    1.00        -14.3541        0.54E-05
         7       2       3    1.00        -12.7122        0.57E-05
         8       2       4    0.00        -10.7393        0.47E-05
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  But, after running the postprocessor, it has been overwritten with:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     15200    ! Version Number
     3   3   2    ! Z, N, 2Mj
  ...
     1      4.1000             711       1    ! Par, WTm, Dim, Npe
  
         1    ! n_states, followed by (i, 2J, nJ, T, -Eb, res)
         1       2       1    0.00        -25.5777        0.65E-05
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  The "basis" files are purely inputs.
  
  The "smwf" files are purely outputs.
  
- `operator application`: These parameters control the application of two-body
  operators.
  
  + `TBMEoperators(1)`: The path to the two-body matrix element (TBME) file.
    Only a single operator may be specified.
    
  + `normalize`: Whether or not to normalize the resultant wave functions.  This
    might be desirable for, e.g., laddering with the J+ operator, to get a new
    normalized state, but undesirable for, e.g., strength function calculations,
    where the norm of the target state contains information on the normalization
    of the strength function.


## HPC run parameters

[Notes on choosing HPC run parameters for mfdn-transitions, as told to mac, by
pjf (07/16/21).]

Q: How many ranks (and threads per rank) should I use?

A: First it is useful to understand how work is split up among ranks.  The bra
and ket are broken up into subvectors, based on the number of "diagonal"
processes in use when they were generated by `mfdn`.  Thus, if there were n
diagonal processes in the `mfdn` run for the ket, and m in the run for the bra,
then the matrix representing the operator in `mfdn-transitions` is effectively
broken into mn tiles, between these pairs of subvectors.  The contributions of
these tiles to the final total matrix element are calculated independently, and
are distributed round-robin over the MPI ranks in `mfdn-transitions`.  While
there are no rigid constraints on choosing the number N of ranks, there are a
few considerations:

- Strong scaling is generally good (except as noted below).

- You want the round-robin distribution to be efficient.  Assuming the
  calculation time to be the same for all tiles, this suggests choosing N as
  a divisor of mn, so that no nodes will be left idle in the final round of
  distribution.  (Corollary: N=mn is the maximum useful number of ranks, as
  N>mn will leave some ranks with no work whatsoever.)

- Calculation time for each tile is typically short (less than a minute)
  compared to the overhead of launching an MPI process (as well as the time to
  run the serial setup codes to generate the OBMEs/TBMEs, if you are doing that
  as part of the same batch job).  So, while choosing N=mn ranks may minimize
  wall time for the actual computation, it will be extremely wasteful of CPU
  time, as most of the billable time will be in MPI setup (and serial setup
  tasks).  For large runs, it is therefore generally better to have N<<mn, to
  minimize MPI overhead relative to actual computation.

- The memory requirements per rank are modest: storage for one ket subvector,
  multiple bra subvectors, and the TBMEs for any two-body operators.  This is
  usually not a limiting factor.

- As usual on HPC systems with a multilayer memory hierarchy, it is probably
  prudent to respect memory boundaries (e.g., sockets or NUMA domains).  Thus,
  as with `mfdn`, you may wish to allocate multiple ranks per node, each
  receiving their proportionate "share" of the total available threads
  efficiently supported by the node (logical cores).

As a rule of thumb, based on the above considerations, N=m or n ranks is a
reasonable choice.

Example: Alice generated Nmax=8 wave functions for 12C, and Bob wants to
calculate transitions.

Alice generated the wave functions using `mfdn` on Cori KNL.  Based on memory
needs, she chose 55 "diagonal" ranks.  That is, she chose to decompose the upper
triangle of the Hamiltonian matrix into 55*(55+1)/2=1540 blocks, and thus ran
mfdn with 1540 ranks. That's all Bob needs to know about Alice's wave function
run, in order to set up his `mfdn-transitions` run.

But let us digress for a moment to review how Alice took the KNL node
configuration into account, as this is useful background for when Bob also takes
the KNL node configuration into account.  Since the KNL nodes are broken into 4
NUMA domains, she split these 1540 ranks over 17 nodes.  With 4x hyperthreading,
and using 64 of the KNL processor's physical cores (to leave a couple in reserve
for system tasks), there are 256 logical cores available per KNL node.  So,
splitting those over 4 ranks, she ran with 64 threads per rank.  Alice was using
the `mcscript` job submission tool `qsubm`, so her submission parameters were
`--ranks=66 --nodes=17 --threads=64` (and `--serialthreads=256` for the "serial",
i.e., OpenMP-parallelized, setup codes from `shell`, which are run before either
`mfdn` or `mfdn-transitions` to generate the OBME/TBME files).

Bob then could therefore reasonably choose anywhere from 1 to 55^2=3025 ranks.
Let's see how this pans out on Cori KNL nodes.  Each node has 64 physical cores
available for the run (to leave a couple in reserve for system tasks), so, with
4x hyperthreading, there are 256 logical cores available per node.

Thus, if Bob ran with one rank per node, he might choose:

   * A single-node run: 1 rank / 1 node / 256 threads per node

   ~~~~
   qsubm ... --ranks=1 --nodes=1 --threads=256 --serialthreads=256
   ~~~~
   
   * Taking as many ranks a diagonal ranks in Alice's run: 55 ranks / 55 nodes /
     256 threads per node

   ~~~~
   qsubm ... --ranks=55 --nodes=55 --threads=256 --serialthreads=256
   ~~~~

However, the memory on each KNL node is not uniformly accessible.  Rather, it is
broken into 4 NUMA domains.  If a rank is spread over multiple NUMA domains,
this may be expected to slow memory access.  So, instead, each rank should be
confined to a single NUMA domain.  It should use only 64 logical cores, which
means we can fit 4 ranks per node.  Thus, Bob might choose:

   * A single-node run: 4 ranks / 1 node / 64 threads per node

   ~~~~
   qsubm ... --ranks=4 --nodes=1 --threads=64 --serialthreads=256
   ~~~~
   
   * Taking as many ranks a diagonal ranks in Alice's run: 55 ranks / 14 nodes /
     64 threads per node

   ~~~~
   qsubm ... --ranks=55 --nodes=14 --threads=64 --serialthreads=256
   ~~~~

   * Or, to again throw 55 nodes at it: 220 ranks / 55 nodes / 64 threads per
     node

   ~~~~
   qsubm ... --ranks=220 --nodes=55 --threads=64 --serialthreads=256
   ~~~~

As an anecdotal example, comparing the timings (runmac0604), we see pretty good
strong scaling from 14 to 55 nodes, if we do it by increasing the number of
ranks (by a factor of 4), but not if we just spread out each rank to use
more threads (and cross over NUMA domains):

   * 55/14/64  => total time with MPI  167.122170000000  => 2338 node-sec
   * 220/55/64 => total time with MPI   50.9975730000000 => 2804 node-sec
   * 55/55/256 => total time with MPI   88.8562030000000 => 4886 node-sec
   * 1/1/256   => total time with MPI 3147.47612600000   => 3147 node-sec


## References

[aktulga2014:mfdn-scalability] "Improving the scalability of symmetric iterative
eigensolver for multi-core platforms", CCPE 26, 2631 (2014).
http://dx.doi.org/10.1002/cpe.3129

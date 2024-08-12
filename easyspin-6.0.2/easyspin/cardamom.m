% cardamom  Trajectory-based simulation of CW-EPR spectra.
%
%   cardamom(Sys,Exp,Par)
%   cardamom(Sys,Exp,Par,Opt)
%   cardamom(Sys,Exp,Par,Opt,MD)
%   spc = cardamom(...)
%   [B,spc] = cardamom(...)
%   [B,spc,TDSignal,t] = cardamom(...)
%   
%   Computes a CW-EPR spectrum of an S=1/2 spin label using stochastic or
%   molecular dynamics trajectories.
%   
%   Sys: stucture with spin system's static and dynamical parameters
%     tcorr          double or numeric vector, size = (1,3)
%                    correlation time (in seconds)
%     logtcorr       double or numeric vector, size = (1,3)
%                    log10 of rotational correlation time (in seconds)
%     Diff           double or numeric vector, size = (1,3)
%                    diffusion rate (rad^2 s^-1)
%     logDiff        double or numeric vector, size = (1,3)
%                    log10 of diffusion rate (rad^2 s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     DiffGlobal     double or numeric vector, size = (1,3)
%                    global diffusion rate (s^-1)
%     Potential      defines an orienting potential using one of the following:
%                    a) set of LMKs and lambdas for an expansion in terms of
%                       Wigner functions [L1 M1 K1 lambda1; L2 M2 K2 lambda2]
%                       lambda can be real or complex
%                    b) 3D array defining the potential over a grid, with the
%                       first dimension representing alpha [0,2*pi], the second
%                       beta [0,pi], and the third gamma [0,2*pi]
%                    c) a function handle for a function that takes three
%                       arguments (alpha, beta, and gamma) and returns the value
%                       of the orientational potential for that orientation. The
%                       function should be vectorized, i.e. work with arrays of
%                       alpha, beta, and gamma.
%     TransRates     numeric, size = (nStates,nStates)
%                    transition rate matrix describing inter-state dynamics
%                    for kinetic Monte Carlo simulations
%     TransProb      numeric, size = (nStates,nStates)
%                    transition probability matrix describing inter-state 
%                    dynamics for kinetic Monte Carlo simulations, note
%                    that a time step must be given to use Sys.TransProb
%                    (alternative input to TransRates; ignored if TransRates
%                    is given)
%     Orientations   numeric, size = (3,nStates)
%                    Euler angles for each state's orientation
%     lw             double or numeric vector, size = (1,2)
%                    vector with FWHM residual broadenings (in mT)
%                         1 element:  GaussianFWHM
%                         2 elements: [GaussianFWHM LorentzianFWHM]
%     lwpp           peak-to-peak linewidths (mT), same format as lw
%                    use either lw or lwpp
%
%   Exp: experimental parameter settings
%     mwFreq         microwave frequency, in GHz (for field sweeps)
%     Range          sweep range, [sweepmin sweepmax], in mT (for field sweep)
%     CenterSweep    sweep range, [center sweep], in mT (for field sweep)
%     nPoints        number of points
%     Harmonic       detection harmonic: 0, 1 (default), 2
%     ModAmp         peak-to-peak modulation amplitude, in mT (field sweeps only)
%     mwPhase        detection phase (0 = absorption, pi/2 = dispersion)
%     Temperature    temperature, in K
%
%   Par: structure with simulation parameters
%     Model      model for spin label dynamics
%                'diffusion': Brownian rotation diffusion with given rotational
%                   diffusion tensor and ordering potential
%                'jump': Markovian jumps between a given set of discrete states
%                'MD-direct': use molecular orientations in MD
%                  trajectories directly as input for simulating the
%                  spectrum
%                'MD-HBD': coarse grain the MD trajectories by using 
%                  the Euler angle probability distribution (for 
%                  pseudopotential) from the spin label's orientations 
%                  to perform further stochastic rotational dynamics 
%                  simulations
%                'MD-HMM': coarse grain the MD trajectories by using 
%                  the spin label's side chain dihedral angles to build 
%                  a hidden Markov model model to perform further 
%                  stochastic jump dynamics simulations
%     dtSpatial  spatial dynamics propagation time step (in seconds)
%                (not used for 'MD-direct')
%     dtSpin     spin dynamics propagation time step (in seconds)
%     nSteps     number of time steps per simulation
%     nTraj      number of trajectories
%     OriStart   numeric, size = (3,1), (1,3), or (3,nTraj)
%                Euler angles for starting orientation(s)
%     nOrients   number of lab-to-molecule orientations to loop over
%     Orients    numeric matrix, size = (nOrients,2)
%                (optional) (phi,theta) angles of lab-to-molecule 
%                orientations. If not given, these are chosen as points
%                on a spherical spiral grid
%
%   Opt: simulation options
%     chkCon         if equal to 1, after the first nSteps of the 
%                    trajectories are calculated, both inter- and intra-
%                    trajectory convergence is checked using the Gelman-
%                    Rubin R statistic such that R<1.1, and if this 
%                    condition is not satisfied, then propagation will be 
%                    extended by either a length of time equal to the 
%                    average of tcorr or by 20% more time steps, whichever 
%                    is larger
%     specCon        if equal to 1, after the first nOrients of the FID
%                    are calculated, both inter- and intra-FID convergence 
%                    are checked using the Gelman-Rubin R statistic such 
%                    that R<1.1, and if this condition is not satisfied, 
%                    then nOrients will be increased by 20% to simulate
%                    additional FIDs until R<1.1 is achieved
%     Verbosity      0: no display, 1: show info
%     Method         string
%                    fast: propagate the density matrix using an 
%                      analytical expression for the matrix exponential in 
%                      the m_S=-1/2 subspace (S=1/2 with up to one nucleus)
%                    ISTOs: propagate the density matrix using
%                      irreducible spherical tensor operators (general, slower)
%     FFTWindow       1: use a Hamming window (default), 0: no window
%     nTrials        number of initialization trials for k-means
%                    clustering; used for the Markov method
%     LagTime        lag time for sliding window processing (only used for
%                    'MD-direct')
%
%   MD: structure with molecular dynamics simulation parameters
%
%     FrameTraj      frame trajectories, size (3,3,nSteps,nTraj)
%     FrameTrajwrtProtein frame trajecoties with respect to protein frame,
%                    size (3,3,nSteps,nTraj)
%     dt             time step (in s) for saving MD trajectory snapshots
%     tLag           time lag (in s) for sampling the MD trajectory to 
%                    determine states and transitions, used for the hidden
%                    Markov model
%     nStates        number of states in the hidden Markov model
%     DiffGlobal     diffusion coefficient for isotropic global rotational
%                    diffusion (s^-1)
%     removeGlobal   integer
%                    1: (default) remove protein global diffusion
%                    0: no removal (e.g. if protein is fixed)
%     LabelName      name of spin label, 'R1' (default) or 'TOAC'
%     HMM            structure, output from 'mdhmm'
%      .TransProb    transition probability matrix
%      .eqDistr      equilibrium distribution vector
%      .mu           center vectors of states
%      .Sigma        covariance matrices of states
%      .viterbiTraj  Viterbi state trajectory (most likely given the dihedrals)
%      .tauRelax     relaxation times of HMM
%      .logLik       log-likelihood of HMM during optimization
%
%   Output:
%     B              magnetic field vector (mT)
%     spc            EPR spectrum
%     TDSignal       FID time-domain signal
%     t              time axis (in s)

% mdhmm    Build hidden Markov model (HMM) from MD trajectory of dihedrals
%
%   HMM = mdhmm(dihedrals,dt,nStates,nLag,Opt)
%
% Input:
%   dihedrals     3D array trajectory of spin-label side chain dihedral angles
%                   (nDihedrals,nTrajectories,nSteps), in radians
%   dt            dihedral trajectory time step, in arbitrary time units
%   nStates       desired number of states for the HMM
%   nLag          desired lag step (number of MD time steps)
%   Opt           structure with options
%     .Verbosity  print to command window if > 0
%     .isSeeded   whether to use rotamer seeds for the centroids in k-means
%     .nTrials    number of trials in k-means clustering (if not seeded)
%
% Output:
%   HMM             structure with HMM parameters
%    .TransProb     transition probability matrix
%    .eqDistr       equilibrium distribution vector
%    .mu            center vectors of states
%    .Sigma         covariance matrices of states
%    .viterbiTraj   Viterbi state trajectory (most likely given the dihedrals)
%    .tauRelax      relaxation times of HMM, in same time units as dt
%    .tLag          lag time, in same time units as dt

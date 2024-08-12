% stochtraj_diffusion  Generate stochastic rotational trajectories
%
%   [t,RTraj] = stochtraj_diffusion(Sys)
%   [t,RTraj,qTraj] = stochtraj_diffusion(...)
%   ... = stochtraj_diffusion(Sys,Par)
%   ... = stochtraj_diffusion(Sys,Par,Opt)
%
%   Sys: stucture with system's dynamical parameters
%
%     tcorr          double or numeric, size = (1,3)
%                    correlation time (in seconds)
%
%     logtcorr       double or numeric, size = (1,3)
%                    log10 of rotational correlation time (in seconds)
%
%     Diff           double or numeric, size = (1,3)
%                    diffusion rate (s^-1)
%
%     logDiff        double or numeric, size = (1,3)
%                    log10 of diffusion rate (s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     Potential      defines an orienting potential using one of the following:
%
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
%
%
%   Par: simulation parameters for Monte Carlo integrator
%
%     nTraj          integer
%                    number of trajectories
%
%     OriStart       numeric, size = (3,1), (1,3), or (3,nTraj)
%                    Euler angles for starting orientation(s) of rotational
%                    diffusion
%
%         When specifying the simulation time provide one of the following
%         combinations:
%         Precedence: t > nSteps,dt > tMax,nSteps
%             
%     t              numeric
%                    array of time points  (in seconds)
%
%     tMax           double
%                    total time of simulation (in seconds)
%
%     dt             double
%                    time step (in seconds)
%
%     nSteps         integer
%                    number of time steps in simulation
%
%
%   Opt: simulation options
%
%     checkConvergence  if set to true, after the first nSteps of the 
%                       trajectories are calculated, both inter- and intra-
%                       trajectory convergence is checked using the Gelman-
%                       Rubin R statistic such that R<1+Opt.convTolerance, and
%                       if this condition is not satisfied, then propagation
%                       will be extended by either a length of time equal to the 
%                       average tcorr or by 20% more time steps, whichever 
%                       is larger.
%
%     convTolerance     Convergence tolerance for Gelman-Rubin R statistic. The
%                       threshold for R is 1 + Opt.convTolerance, e.g. if 
%                       Opt.convTolerance = 1e-6, then the threshold R is
%                       1.000001.
%
%     Verbosity         0: no display, 1: show info
%
%
%   Output:
%
%     t              matrix, size = (nSteps,1) 
%                    time points of the trajectory (in seconds)
%
%     RTraj          4D array, size = (3,3,nSteps,nTraj)
%                    trajectories of rotation matrices
%
%     qTraj          3D array, size = (4,nSteps,nTraj)
%                    trajectories of normalized quaternions

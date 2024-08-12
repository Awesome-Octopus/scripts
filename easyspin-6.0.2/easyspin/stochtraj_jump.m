% stochtraj_jump  Generate stochastic trajectories of Markovian jumps using
%                 kinetic Monte Carlo
%
%   [t,RTraj] = stochtraj_jump(Sys)
%   ... = stochtraj_jump(Sys,Par)
%   ... = stochtraj_jump(Sys,Par,Opt)
%   [t,RTraj,qTraj] = stochtraj_jump(...)
%   [t,RTraj,qTraj,stateTraj] = stochtraj_jump(...)
%   [t,stateTraj] = stochtraj_jump(...)
%
%   Sys: stucture with system's dynamical parameters
%
%     TransRates     numeric, size = (nStates,nStates)
%                    transition rate matrix describing inter-state dynamics
%                    for kinetic Monte Carlo simulations
%
%     TransProb      numeric, size = (nStates,nStates)
%                    transition probability matrix describing inter-state 
%                    dynamics for kinetic Monte Carlo simulations, note
%                    that a time step must be given to use Sys.TransProb
%                    (alternative input to TransRates; ignored if TransRates
%                    is given)
%
%     Orientations   numeric, size = (3,nStates)
%                    Euler angles for each state's orientation
%
%
%   Par: simulation parameters for Monte Carlo integrator
%
%     nTraj          integer
%                    number of trajectories
%
%     StatesStart    numeric, size = (1,1) or (nTraj,1)
%                    starting states of trajectories
%
%         When specifying the simulation time provide one of the following
%         combinations:
%         Precedence: nSteps,dt > tMax,nSteps
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
%     statesOnly        1 or 0
%                       specify whether or not to use a discrete model purely
%                       to generate states and not quaternion orientations
%
% %     checkConvergence  if equal to 1, after the first nSteps of the 
% %                       trajectories are calculated, both inter- and intra-
% %                       trajectory convergence is checked using the Gelman-
% %                       Rubin R statistic such that R<1.1, and if this 
% %                       condition is not satisfied, then propagation will be 
% %                       extended by either a length of time equal to the 
% %                       average of tcorr or by 20% more time steps, whichever 
% %                       is larger
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
%
%     stateTraj      2D array, size = (nSteps,nTraj)
%                    trajectories of states

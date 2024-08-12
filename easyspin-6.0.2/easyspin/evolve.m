% evolve  Time-domain evolution of density matrix
%
%   td = evolve(Sig,Det,Ham,n,dt);
%   td = evolve(Sig,Det,Ham,n,dt,IncScheme);
%   td = evolve(Sig,Det,Ham,n,dt,IncScheme,Mix);
%
%   Evolves the density matrix Sig under the Hamiltonian Ham with time
%   step dt n-1 times and detects using Det after each step. Hermitian
%   input matrices are assumed. Ham is assumed to be in units of MHz.
%
%   n gives the number of points along each dimension (a number for 1D
%   experiments, a 2-vector for 2D experiments). dt contains the time increment
%   along each dimension, in microseconds.
%
%   IncScheme determines the incrementation scheme and can be one of the
%   following (up to five incrementation periods, up to two dimensions)
%
%     1D experiments:
%       [1]           simple FID, 3p-ESEEM, echo transient, DEFENCE
%       [1 1]         2p-ESEEM, CP, 3p- and 4p-RIDME
%       [1 -1]        3p-DEER, 4p-DEER, PEANUT, 5p-RIDME
%       [1 -1 1 -1]
%       [1 1 -1 -1]   SIFTER
%       [1 -1 -1 1]   7p-DEER
%
%     2D experiments:
%       [1 2]         3p-ESEEM echo transient, HYSCORE, DONUT-HYSCORE
%       [1 1 2]       2p-ESEEM etc. with echo transient
%       [1 -1 2]      3p-DEER, 4p-DEER etc. with echo transient
%       [1 2 1]       2D 3p-ESEEM
%       [1 1 2 2]     2D refocused 2p-echo
%       [1 2 2 1]     2D CP
%       [1 2 -2 1]    2D PEANUT
%       [1 1 -1 -1 2] SIFTER with echo transient
%       [1 -1 -1 1 2] 7p-DEER with echo transient
%
%   [1] is the default. 1 and -1 indicate incrementation and decrementation
%   along the first dimensions, and 2 and -2 are analogous for the second
%   dimensions. For more details, see the documentation.
%
%   Mix is a cell array containing the propagators of the mixing
%   block(s). It is required for experiments with more than 1 sweep period.
%
%   td is the output signal, a vector for 1D experiments and a matrix for
%   2D experiments. td(1) is the value obtained by detecting Sig without evolution.

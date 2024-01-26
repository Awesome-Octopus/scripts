% exciteprofile      Excitation profile calculation for arbitrary pulses
%
% [offsets,Mag] = exciteprofile(t,IQ)
% [offsets,Mag] = exciteprofile(t,IQ,offsets)
%
% Input:
%   t           = time axis for defined waveform in �s
%   IQ          = in-phase and quadrature part of the pulse function
%   offsets     = axis of frequency offsets in MHz for which to compute the
%                 excitation profile (default: approximately �2*BW centered
%                 at the pulse center frequency, 201 points)
% Output:
%   offsets     = axis of frequency offsets in MHz
%   Mag         = excitation profile as Mi/M0, i = x,y,z in array
%

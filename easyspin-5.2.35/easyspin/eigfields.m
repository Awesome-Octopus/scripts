% eigfields  Resonance fields by the eigenfield method 
%
%   B = eigfields(Sys, Par,)
%   B = eigfields(Sys, Par, Opt)
%   [B, Int] = eigfields(...)
%
%   Calculates all resonance fields of spin system
%   Sys solving a generalized eigenvalue problem
%   in Liouville space together with transition pro-
%   babilitites.
%
%   Input:
%   - Sys: spin system specification structure
%   - Par: structure with fields
%        mwFreq - spectrometer frequency [GHz]
%        Mode - 'parallel' or 'perpendicular' (default)
%          direction of mirowave field relative to static field
%        Range - [Bmin Bmax] If set, compute only eigenfields
%           between Bmin and Bmax. [mT]
%   - Opt: options structure with fields
%        Threshold - if set, return only transitions with
%          relative intensity above Threshold.
%
%   Output:
%   - B:   cell array of all resonance fields [mT]
%   - Int: transition intensities [MHz^2/mT^2]

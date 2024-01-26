% nucspinkeep  Remove nuclear spins from spin system
%
%   NewSys = nucspinrmv(Sys,keepidx)
%
%   Removes one or more nuclei from the spin system
%   Sys and returns the result in NewSys. keepidx a vector
%   of nuclear indices as they appear in Sys. E.g. keepidx=1
%   removes all except the first nucleus, and keepidx=[2 4]
%   removes all except the second and the fourth nucleus.
%
%   Example:
%    Sys.Nucs = '13C,1H,14N,14N';
%    Sys.A = [14 5 7 8.2];
%    Sys = nucspinkeep(Sys,2);  % keeps only 1H

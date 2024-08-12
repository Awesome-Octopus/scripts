% nucspinrmv   Remove nuclear spin from spin system
%
%   NewSys = nucspinrmv(Sys,rmvidx)
%
%   Removes one or more nuclei from the spin system
%   Sys and returns the result in NewSys. rmvidx a vector
%   of nuclear numbers as they appear in Sys. E.g. rmvidx=1
%   removes the first nucleus, rmvidx=[2 4] removes the
%   second and the fourth.
%
%   Example:
%    Sys = struct('S',1/2,'g',[2 2 3],'Nucs','14N,1H','A',[4 6 10; 1 1 2]);
%    Sys = nucspinrmv(Sys,1);  % removes the 1H

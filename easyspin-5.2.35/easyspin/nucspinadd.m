% nucspinadd  Adds a nuclear spin to a spin system
%
%    NewSys = nucspinadd(Sys,Nuc,A)
%    NewSys = nucspinadd(Sys,Nuc,A,AFrame)
%    NewSys = nucspinadd(Sys,Nuc,A,AFrame,Q)
%    NewSys = nucspinadd(Sys,Nuc,A,AFrame,Q,QFrame)
%
%    NewSys = nucspinadd(Sys,Nuc,Afull)
%    NewSys = nucspinadd(Sys,Nuc,Afull,Qfull)
%
%    Add the nuclear isotope Nuc (e.g. '14N') to
%    the spin system structure, with the hyperfine
%    values A, the hyperfine tilt angles AFrame, the
%    quadrupole values Q and the quadrupole tilt angles
%    QFrame. Any missing parameter is assumed to be [0 0 0].
%
%    Alternatively, full 3x3 hyperfine and quadrupole
%    matrices can be specified in Afull and Qfull.
%
%    Examples:
%     Sys = struct('S',1/2,'g',[2 2 2.2]);
%     Sys = nucspinadd(Sys,'Cu',[50 50 520]);
%     Sys = nucspinadd(Sys,'14N',[20 0 0; 0 30 0; 0 0 50]);

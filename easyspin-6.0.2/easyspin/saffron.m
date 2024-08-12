% saffron    Simulate pulse EPR signals and spectra
%
%     S = saffron(Sys,Exp,Opt)
%     [x,S] = saffron(Sys,Exp,Opt)
%     [x,S,info] = saffron(Sys,Exp,Opt)
%
%  Inputs:
%     Sys   ... spin system with electron spin and ESEEM nuclei
%     Exp   ... experimental parameters (time unit µs)
%     Opt   ... simulation options
%
%  Outputs:
%     x       ... axis/axes for S, contains all indirect dimensions and for
%                 transient detection transient time axis of transient, for
%                 ENDOR the frequency axis
%     S       ... simulated signal (ESEEM) or spectrum (ENDOR)
%     info    ... structure with FFT of ESEEM signal etc.

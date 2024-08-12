% spidyan    Simulate spindyanmics during pulse EPR experiments
%
%     [t,Signal] = spidyan(Sys,Exp,Opt)
%     [t,Signal,info] = spidyan(Sys,Exp,Opt)
%
%     Inputs:
%       Sys   ... spin system with electron spin and ESEEM nuclei
%       Exp   ... experimental parameters (time unit µs)
%       Opt   ... simulation options
%
%     Outputs:
%       t                 ... time axis
%       Signal            ... simulated signals of detected events
%       info              ... structure with fields:
%         .FinalState        ... density matrix/matrices at the end of the
%                                experiment
%         .StateTrajectories ... cell array with density matrices from each
%                                timestep during evolution
%         .Events            ... structure containing events

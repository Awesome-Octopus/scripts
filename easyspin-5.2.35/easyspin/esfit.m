% esfit   least-squares fitting for EPR spectral simulations
%
%   esfit(simfunc,expspc,Sys0,Vary,Exp)
%   esfit(simfunc,expspc,Sys0,Vary,Exp,SimOpt)
%   esfit(simfunc,expspc,Sys0,Vary,Exp,SimOpt,FitOpt)
%   bestsys = esfit(...)
%   [bestsys,bestspc] = esfit(...)
%
%     simfunc     simulation function name ('pepper', 'garlic', 'salt', ...
%                   'chili', or custom function), or function handle
%     expspc      experimental spectrum, a vector of data points
%     Sys0        starting values for spin system parameters
%     Vary        allowed variation of parameters
%     Exp         experimental parameter, for simulation function
%     SimOpt      options for the simulation algorithms
%     FitOpt      options for the fitting algorithms
%        Method   string containing kewords for
%          -algorithm: 'simplex','levmar','montecarlo','genetic','grid',
%                      'swarm'
%          -target function: 'fcn', 'int', 'dint', 'diff', 'fft'
%        Scaling  string with scaling method keyword
%          'maxabs' (default), 'minmax', 'lsq', 'lsq0','lsq1','lsq2','none'
%        OutArg   two numbers [nOut iOut], where nOut is the number of
%                 outputs of the simulation function and iOut is the index
%                 of the output argument to use for fitting

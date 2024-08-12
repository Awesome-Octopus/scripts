% esfit_montecarlo   Monte Carlo algorithm for least-squares fitting
%
%    x = esfit_montecarlo(fcn,lb,ub,FitOpt)
%
%    fcn ...  scalar function to minimize, f(x), where x is an array
%    lb   ... lower bounds for x
%    ub   ... upper bounds for x
%    FitOpt ... options
%       nTrials       number of trials
%       maxTime       maximum time to run, in minutes
%       Verbosity    1, if progress information should be printed
%       TolFun        error threshold below which fitting stops

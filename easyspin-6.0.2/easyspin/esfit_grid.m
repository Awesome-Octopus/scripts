% esfit_grid   Function minimization using grid search
%
%  xfit = esfit_grid(fcn,lb,ub)
%  ...  = esfit_grid(fcn,lb,ub,FitOpt)
%  [xfit,info] = ...
%
%  Finds x that minimzes fcn(x), running over a grid of parameter values.
%
%  Input:
%    fcn     function handle of fcn(x) to minimize
%    lb      lower bounds of parameters
%    ub      lower bounds of parameters
%    FitOpt  structure with algorithm parameters
%      .TolFun
%      .GridSize
%      .randGrid
%      .maxGridPoints
%      .IterFcn   function that is called after each iteration

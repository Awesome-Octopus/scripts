%esfit_swarm    
%
%  xfit = esfit_swarm(fcn,lb,ub)
%  ...  = esfit_swarm(fcn,lb,ub,FitOpt)
%  [xfit,info] = ...
%
%  ...
%
%  Input:
%    fcn     function handle of fcn(x) to minimize
%    lb      lower bounds of parameters
%    ub      lower bounds of parameters
%    Opt     structure with algorithm parameters
%      .IterFcn   function that is called after each iteration

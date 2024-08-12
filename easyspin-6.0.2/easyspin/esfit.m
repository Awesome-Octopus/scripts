% esfit   Least-squares fitting for EPR and other data
%
%   esfit(data,fcn,p0,vary)
%   esfit(data,fcn,p0,lb,ub)
%   esfit(___,FitOpt)
%
%   pfit = esfit(___)
%
% Input:
%     data        experimental data, a vector of data points or a cell
%                   array of datasets for global fitting 
%     fcn         simulation/model function handle (@pepper, @garlic, ...
%                   @salt, @chili, or handle to user-defined function)
%                   a user-defined fcn should take a parameter vector p
%                   and return simulated data datasim: datasim = fcn(p)
%     p0          starting values for parameters
%                   EasySpin-style functions: {Sys0,Exp0} or {Sys0,Exp0,Opt0}
%                   other functions: n-element vector
%     vary        allowed variation of parameters
%                   EasySpin-style functions: {vSys} or {vSys,vExp} or {vSys,vExp,vOpt}
%                   other functions: n-element vector
%     lb          lower bounds of parameters
%                   EasySpin-style functions: {lbSys,lbExp} or {lbSys,lbExp,lbOpt}
%                   other functions: n-element vector
%     ub          upper bounds of parameters
%                   EasySpin-style functions: {ubSys,ubExp} or {ubSys,ubExp,ubOpt}
%                   other functions: n-element vector
%     FitOpt      options for esfit
%        .Method  string containing keywords for
%           -algorithm: 'simplex','levmar','montecarlo','genetic','grid','swarm'
%           -target function: 'fcn', 'int', 'dint', 'diff', 'fft'
%        .AutoScale 'lsq', 'maxabs', 'none'; default 'lsq' for
%                 EasySpin simulation functions, otherwise 'none'
%        .BaseLine 0, 1, 2, 3 or [] (or vector for global fitting with different 
%                 baseline order for different datasets)
%        .OutArg  two numbers [nOut iOut], where nOut is the number of
%                 outputs of the simulation function and iOut is the index
%                 of the output argument to use for fitting
%        .Mask    array of 1 and 0 the same size as data vector
%                 values with mask 0 are excluded from the fit 
%                 (cell array for data input consisting of multiple datasets)
%        .weight  array of weights to use when combining residual vectors
%                 of all datasets for global fitting
% Output:
%     fit           structure with fitting results
%       .pfit       fitted parameter vector (contains only active fitting parameters)
%       .pnames     variable names of the fitted parameters
%       .pfit_full  parameter vector including inactive fitting parameters (in GUI)
%       .argsfit    fitted input arguments (if EasySpin-style)
%       .pstd       standard deviation for all parameters
%       .ci95       95% confidence intervals for all parameters
%       .cov        covariance matrix for all parameters
%       .corr       correlation matrix for all parameters
%       .p_start    starting parameter vector for fit
%       .fitraw     fit, as returned by the simulation/model function
%       .fit        fit, including the fitted scale factor
%       .scale      fitted scale factor
%       .baseline   fitted baseline
%       .mask       mask used for fitting
%       .residuals  residuals
%       .ssr        sum of squared residuals
%       .rmsd       root-mean square deviation between input data and fit
%       .bestfithistory  structure containing a list of fitting parameters
%                        corresponding to progressively improved rmsd
%                        values during fitting process and corresponding
%                        rmsd values, for EasySpin functions, a conversion
%                        function returning the EasySpin input structures
%                        given a selected set of fitting parameters is also
%                        included
%

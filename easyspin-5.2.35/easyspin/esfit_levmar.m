%esfit_levmar  Levenberg-Marquardt nonlinear least squares fitting
%
%   xfit = esfit_levmar(funfcn,x0)
%   [xfit,Info] = esfit_levmar(funfcn,x0)
%   ... = esfit_levmar(funfcn,x0,Opt)
%   ... = esfit_levmar(funfcn,x0,Opt,p1,p2,...)
%
%   Find  xm = argmin{F(x)} , where  x = [x_1, ..., x_n]  and
%   F(x) = sum(f_i(x)^2)/2. The functions  f_i(x) (i=1,...,m)
%   must be given by a Matlab function with declaration
%              function  f = funfcn(x,p1,p2,...)
%   p1,p2,... are parameters of the function and can be of any type and size.
%
%   The parameter search range is restricted to -1...+1 along each
%   dimension.
%
% Input
%   funfcn Handle to the function.
%   x0     Starting vector in parameter space
%   Opt    Options structure
%            lambda    starting value of Marquardt parameter
%            Gradient  termination threshold for gradient
%            TolStep   termination threshold for step
%            maxTime   termination threshold for time
%            delta     step width for Jacobian approximation
%   p1,p2,... are passed directly to the function funfcn.    
%
% Output
%   xfit    Converged vector in parameter space
%   Info    Performance information, vector with 7 elements:
%           info(1) = final value of F(x)
%           info(2) = final value of ||F'||inf
%           info(3) = final value of ||dx||2
%           info(4) = final value of mu/max(A(i,i)) with A = Je'* Je
%           info(5) = number of iterations
%           info(6) = 1 :  stopped by small gradient
%                     2 :  stopped by small x-step
%                     3 :  max no. of iterations exceeded 
%                    -4 :  dimension mismatch in x, f, B0
%                    -5 :  overflow during computation
%                    -6 :  error in approximate Jacobian
%           info(7) = number of function evaluations

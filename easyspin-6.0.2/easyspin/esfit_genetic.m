% esfit_genetic   Genetic algorithm for least-squares fitting
%
%    x = esfit_genetic(fcn,nParams,FitOpt)
%
%    fcn      ... scalar function to minimize
%    nParams  ... number of parameters
%    FitOpt   ... options
%       .PopulationSize   number of individuals per generation
%       .EliteCount       number of elite individuals
%       .maxGenerations   maximum number of generations
%       .Verbosity       1, if progress information should be printed
%       .TolFun           error threshold below which fitting stops

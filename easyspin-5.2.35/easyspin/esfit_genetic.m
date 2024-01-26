% esfit_genetic   Genetic algorithm for least-squares fitting
%
%    x = esfit_genetic(funfcn,nParams,FitOpt,varargin)
%
%    funfcn ... scalar function to minimize
%    nParams ... number of parameters
%    FitOpt ... options
%       PopulationSize   number of individuals per generation
%       EliteCount       number of elite individuals
%       Range            parameter range (from -Range to +Range)
%       maxGenerations   maximum number of generations
%       PrintLevel          1, if progress information should be printed
%       TolFun           error threshold below which fitting stops

% Fitting a two-component spectrum
%====================================================================
clear, clc

% Since we don't have any experimental data available, we create a
% two-component spectrum and add noise. If you use this example as
% a basis for your fittings, replace this code by code that loads
% your experimental data.

Sys1.g = [2.00703 2.00308, 2.00131];
%Sys1.lwpp = 0.1;  % mT
Sys1.Nucs = '14N';
Sys1.A = [17.7, 83.5];
Sys1.weight = 1;
Sys1.logtcorr = -8;
%Sys1.lw = 0.1;
% Sys2.g = [2.008 2.006 2.003];
% Sys2.lwpp = 2;  % mT
% Sys2.weight = 0.3;
% Sys2.logtcorr = -9

[B,spc, params] = eprload("~/EMX/sl_lipids/062424_12_doxyl_DPPC_307K.DSC");
B = B/10;
% spc = addnoise(spc,50,'n');
% baseline subtraction
%spc = spc - mean(spc(1:100))

Exp.mwFreq = params.MWFQ/1e9;  % GHz
Exp.Range = [324.28, 335.28];  % mT

% Next we set up the least-squares fitting.
% First comes a starting set of parameters (which we obtain by reusing
% the spin systems from the simulation and changing a few values)
Sys = {Sys1};

% Next, we specify which parameter we want to be fitted and by how much
% the fitting algorithm can vary it.
Vary1.g = [0.001 0.001 0.001];
Vary1.logtcorr = 4;
Vary1.A = [40 40];
%Vary1.lw = 0.05;
%Vary1.lwpp = 0.02;
%Vary1.weight = 500
%Vary2.lwpp = 0.9;
%Vary2.logtcorr = 9
%Vary2.weight = 500;
Vary = {Vary1};


% Call the fitting function
SimOpt.Method = 'perturb';
SimOpt.AutoScale = 'lsq';
FitOpt.Method = 'simplex'; % simplex algorithm, integrals of spectra
FitOpt.maxTime = 1;
FitOpt.TolEdgeLength = 1e-3;
esfit_legacy(spc,@chili,{Sys,Exp,SimOpt},{Vary},FitOpt);



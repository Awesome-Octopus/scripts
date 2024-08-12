% Fitting a two-component spectrum
%====================================================================
% clear
% 
% % Since we don't have any experimental data available, we create a
% % two-component spectrum and add noise. If you use this example as
% % a basis for your fittings, replace this code by code that loads
% % your experimental data.
% 
% sys1.g = [2.00468 2.01145 2.00343];
% sys1.A = [13.5 13.5 99];
% 
% sys1.logtcorr = -8.29035;
% sys1.Nucs = '14N';
% sys1.weight = 114.763;
% sys1.lwpp = 0.0059870;
% 
% sys2.g = [2.00468 2.01145 2.00343];
% sys2.A = [13.5 13.5 99];
% 
% sys2.logtcorr = -8.29035;
% sys2.Nucs = '14N';
% sys2.weight = 114.763;
% sys2.lwpp = 0.0059870;
% 
% 
% 
% [B,spc, params] = eprload();
% exp.mwFreq = params.MWFQ/(1e9);
% B = B/10;
% exp.Range = [min(B), max(B)];

% % when we've found a good fit save it, then comment out the stuff above and
% % start from the saved parameters by uncommenting the stuff below.
% 
Sys1 = fit1.Sys{1};
Sys2 = fit1.Sys{2};


% Next we set up the least-squares fitting.
% First comes a starting set of parameters (which we obtain by reusing
% the spin systems from the simulation and changing a few values)



% Next, we specify which parameter we want to be fitted and by how much
% the fitting algorithm can vary it.


% vary1.weight = 150;
% vary1.logtcorr = 1.5;
% vary1.g = [0.015 0.03 0.015];
% vary1.lwpp = 0.01;
% vary1.A = [2 2 5];
% 
% vary2.g = [0.015 0.03 0.015];
% vary2.A = [2 2 5];
% vary2.lwpp = 0.01;
% vary2.logtcorr = 1.5;
% vary2.weight = 150;


% Calling the fitting function
% SimOpt.Method = 'perturb';
% FitOpt.Method = 'simplex int'; % simplex algorithm, integrals of spectra
% 
% 
% esfit('chili',spc,{sys1,sys2},{vary1,vary2},exp,SimOpt,FitOpt);

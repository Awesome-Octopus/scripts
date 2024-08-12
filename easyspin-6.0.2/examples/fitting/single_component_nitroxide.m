
%====================================================================
% fits simulation parameters for a nitroxide with one motional component
% 
function fit = single_component_nitroxide (B, spc, params, starting_fit)

    Sys.Nucs = '14N';
    B = B/10; %experimental field is given in gauss, but needs to be millitesla
    Exp.mwFreq = params.MWFQ/1e9;  % GHz
    Exp.Range = [params.XMIN/10, (params.XMIN/10 + params.XWID/10)];  % mT
    extended_fit = false;
    % start from the fit of the previous spectrum in the sequence
    % only fit the logtcorr
    if ~isempty(starting_fit)

        Sys.g = starting_fit.argsfit{1,1}.g;
        Sys.logtcorr = starting_fit.argsfit{1,1}.logtcorr;
        Sys.A = starting_fit.argsfit{1,1}.A;
        FitOpt.maxTime = 1;

        % if you tell it to it will adjust g and A values each time. this
        % makes the whole process more accurate but takes considerably longer.
        if extended_fit == true
            Vary.A = [2, 2, 4];
            %Vary.g = [0.0005, 0.001, 0.0005];
            
            FitOpt.maxTime = 5;
        end
        FitOpt.Method = 'simplex';     
        Vary.logtcorr = 0.3;
        Vary.A = [2, 2, 6];

    % the intial guess is blind and many parameters are tweaked
    % starting from reasonable values for nitroxides
    else

        Sys.g = [2.0086, 2.0081, 2.0033];
        Sys.A = [23.6717, 18.6602, 74.5816];
        Sys.logtcorr = [-7.93959, -5.52869];
        Vary.g = [0.001, 0.001, 0.001];
        Vary.A = [3, 3, 10];
        FitOpt.Method = 'montecarlo';
        FitOpt.nTrials = 1000;
        Vary.logtcorr = [1, 1];
   

        % Call the fitting function
        SimOpt.Method = 'perturb';
        SimOpt.AutoScale = 'lsq';
        fit = esfit_legacy(spc,@chili,{Sys,Exp,SimOpt},{Vary},FitOpt);
        Sys.g = fit.argsfit{1,1}.g;
        Sys.logtcorr = fit.argsfit{1,1}.logtcorr;
        Sys.A = fit.argsfit{1,1}.A;
        Vary.g = [0.001 0.001 0.001];
        Vary.logtcorr = [2, 2];
        Vary.A = [5 5 15];
        FitOpt.Method = 'simplex';
        FitOpt.baseline = 'baseline';
        FitOpt.maxTime = 5;

    end

    

    % Call the fitting function
    SimOpt.Method = 'perturb';
    SimOpt.AutoScale = 'lsq';
    
    fit = esfit_legacy(spc,@chili,{Sys,Exp,SimOpt},{Vary},FitOpt);

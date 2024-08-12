
%====================================================================
% fits simulation parameters for a nitroxide with one motional component
% 
function fit = two_component_nitroxide (varargin)

    Sys1.Nucs = '14N';
    Sys2.Nucs = '14N';
    try
        [B, spc, params, fits] = deal(varargin{:});
        B = B/10; %experimental field is given in gauss, but needs to be millitesla
        Exp.mwFreq = params.MWFQ/1e9;  % GHz
        Exp.Range = [params.XMIN/10, (params.XMIN/10 + params.XWID/10)];  % mT

        if length(varargin == 3) | isempty(fits)
            error('no existing fit supplied, starting a priori')
            
        end

        % start from the fit of the previous spectrum in the sequence

        Sys1.g = fit1.argsfit{1,1}.g;
        Sys1.logtcorr = fit1.argsfit{1,1}.logtcorr;
        Sys1.A = fit1.argsfit{1,1}.A;
        % Vary1.g = [0.001 0.002 0.001];
        % Vary1.A = [5 5 10];
        Vary1.logtcorr = 2;

        Sys2.g = fit1.argsfit{2,1}.g;
        Sys2.logtcorr = fit1.argsfit{2,1}.logtcorr;
        Sys2.A = fit1.argsfit{2,1}.A;
        % Vary2.g = [0.001 0.002 0.001];
        % Vary2.A = [5 5 10];
        Vary2.logtcorr = 2;
  

        FitOpt.Method = 'simplex';     
        FitOpt.maxTime = 5;

    catch ME

        try % when the fit is not provided

            % the intial guess is blind and many parameters are tweaked
            % starting from reasonable values for nitroxides

            Sys1.g = [2.008 2.008, 2.004];
            Sys1.A = [18, 18, 86];
            Sys1.logtcorr = -9;
            Sys1.weight = 0.5;
    
            Vary1.g = [0.001 0.002 0.001];
            Vary1.A = [10 10 20];
            Vary1.weight = 0.5;
            Vary1.logtcorr = 2;
    
            Sys2.g = [2.009 2.009, 2.005];
            Sys2.A = [18, 18, 86];
            Sys2.logtcorr = -9;
            Sys2.weight = 0.5;
    
            Vary2.g = [0.001 0.002 0.001];
            Vary2.A = [10 10 20];
            Vary2.weight = 0.5;
            Vary2.logtcorr = 2;
    
    
            FitOpt.Method = 'montecarlo';
            FitOpt.nTrials = 500;
    
            Sys = {Sys1, Sys2};
            Vary = {Vary1, Vary2};
    
            % Call the fitting function
            SimOpt.Method = 'perturb';
            SimOpt.AutoScale = 'lsq';

            if length(varargin) == 0 | isinteger(varargin{1}) 
                error('No input data structure given, operating in interactive mode')
            end


            fit = esfit_legacy(spc,@chili,{Sys,Exp,SimOpt},{Vary},FitOpt);
    
            Sys1.g = fit.argsfit{1,1}.g;
            Sys1.logtcorr = fit.argsfit{1,1}.logtcorr;
            Sys1.A = fit.argsfit{1,1}.A;
            Vary.g = [0.001 0.002 0.001];
            Vary.logtcorr = 2;
            Vary.A = [10 10 20];
    
            Sys2.g = fit.argsfit{2,1}.g;
            Sys2.logtcorr = fit.argsfit{2,1}.logtcorr;
            Sys2.A = fit.argsfit{2,1}.A;
            Vary2.g = [0.001 0.002 0.001];
            Vary2.logtcorr = 2;
            Vary2.A = [10 10 20];
    
            FitOpt.Method = 'simplex';
            FitOpt.maxTime = 10;
    
        
            SimOpt.Method = 'perturb';
            SimOpt.AutoScale = 'lsq';
            

            
            fit = esfit_legacy(spc,@chili,{{Sys1, Sys2},Exp,SimOpt},...
                {Vary1, Vary2},FitOpt);


        catch ME % if no input spectrum is provided
            [temp_and_field, data, params] = eprload();
            field = temp_and_field{1};
            temps = temp_and_field{2};
            B = field/10; %experimental field is given in gauss, but needs to be millitesla
            Exp.mwFreq = params.MWFQ/1e9;  % GHz
            Exp.Range = [params.XMIN/10, (params.XMIN/10 + params.XWID/10)];  % mT

            if isempty(varargin)
                
                v = floor(numel(temps)/2);
                spc = data(:,v);

            elseif isinteger(varargin{1})
                spc = data(:,varargin{1});
            end
                
        
            % Call the fitting function
            SimOpt.Method = 'perturb';
            SimOpt.AutoScale = 'lsq';
            esfit_legacy(spc,@chili,{Sys,Exp,SimOpt},{Vary},FitOpt);
        end
    end

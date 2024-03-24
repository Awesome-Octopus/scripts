%% Powersat Multifit
function condition = powersat_multifit

condition = powersat_multiload;


    %% set this value to 1 to allow masking of low attenuation points from data analysis
    doDataMask = 0;


    %% mask bad data points if that was flagged
    if doDataMask == 1

        numMaskPts = input("Mask the highest X power steps: X = ");

        for p = 1:numel(condition)

            %then for each run ...
            for i =1:numel(condition(p).rep)

                %find how many decibels you are stepping by based on your
                %minimum power and number of (unmasked) power steps
                condition(p).rep(i).attstep = (10*log10(200/condition(p).rep(i).params.YMIN)...
                + 10*log10(200/(condition(p).rep(i).params.YMIN + ...
                condition(p).rep(i).params.YWID)))...
                /condition(p).rep(i).params.YPTS;


                powermaxatt = 10*log10(200/(condition(p).rep(i).params.YWID - ...
                condition(p).rep(i).params.YMIN));

                % change the YPTS in params to the number of attenuation steps -
                % masked steps
                condition(p).rep(i).params.YPTS = condition(p).rep(i).params.YPTS - numMaskPts;

                % and delete the first X columns in the y data
                condition(p).rep(i).ydata = condition(p).rep(i).ydata(:,numMaskPts+1:end);

                % and delete the first X rows in xdata
                condition(p).rep(i).xdata{1,2} = condition(p).rep(i).xdata{1,2}(numMaskPts+1:end);

                % scale the YWID param to adjust max power point
                condition(p).rep(i).params.YWID = 200/10^((powermaxatt+numMaskPts*...
                condition(p).rep(i).attstep)/10)-condition(p).rep(i).params.YMIN;

            end    
        end    


    end
    %%
    for i = 1:numel(condition) % for each condition
        for p = 1:numel(condition(i).rep) % for each replicate

            [m, n] = size(condition(i).rep(p).ydata);

            condition(i).rep(p).ysmooth = datasmooth(condition(i).rep(p).ydata,10);
            % smooth the y data

            for q = 1:condition(i).rep(p).params.YPTS 
            % and for each power step within each replicate within each condition


                [condition(i).rep(p).maxval(q), condition(i).rep(p).maxpnt(q)] = ...
                    max(condition(i).rep(p).ysmooth(:,q));
                [condition(i).rep(p).minval(q), condition(i).rep(p).minpnt(q)] =... 
                    min(condition(i).rep(p).ysmooth(:,q)); 
                % detemine the peak to trough height 
                % for each scan of the smoothed data

                condition(i).rep(p).height(q) = condition(i).rep(p).maxval(q)- ...
                    condition(i).rep(p).minval(q);

                condition(i).rep(p).delH(q) = condition(i).rep(p).xdata{1,1}(condition(i).rep(p).minpnt(q)) ...
                    - condition(i).rep(p).xdata{1,1}(condition(i).rep(p).maxpnt(q));




            end


            if condition(i).rep(p).height(1) > condition(i).rep(p).height(n)
                 condition(i).rep(p).height = flip(condition(i).rep(p).height);
            end


            %%Extract the sqrt of the power each attentuation step was taken at
            powerminatt = 10*log10(200/(condition(i).rep(p).params.YMIN));
            powermaxatt = 10*log10(200/(condition(i).rep(p).params.YWID - ...
                condition(i).rep(p).params.YMIN));
            poweratt = linspace(powerminatt,powermaxatt,n);
            condition(i).rep(p).power = 200./(10.^(poweratt/10));
            condition(i).rep(p).sqrtpower = condition(i).rep(p).power.^0.5;



            %% Fit the Data to obtain P1/2
            options = optimset('MaxFunEvals',10000'); 
            % Sets the Maximum number of function evaluations

            [condition(i).rep(p).x,condition(i).rep(p).fvalA] = ...
                fmincon(@powerfit,[max(condition(i).rep(p).height)...
                /4,10,1.2],[],[],[],[],[0,0,0.5],[1e6,1e6,1.5],[],options,...
                condition(i).rep(p).power, condition(i).rep(p).sqrtpower,...
                condition(i).rep(p).height); 
             % Uses the Matlab builtin fminsearch to minimize the t1fit function

            [condition(i).rep(p).y,condition(i).rep(p).fvalB] = ...
                fminsearch(@linearfit,[(condition(i).rep(p).height(4)-...
                condition(i).rep(p).height(1))/(condition(i).rep(p).sqrtpower(4)...
                -condition(i).rep(p).sqrtpower(1)), condition(i).rep(p).height(1)],...
                options, condition(i).rep(p).sqrtpower, condition(i).rep(p).height);

            condition(i).rep(p).xdata = linspace(min(condition(i).rep(p).sqrtpower),...
                max(condition(i).rep(p).sqrtpower),1001); 
            % creates x data array for the plotted fit

            [m n] = size(condition(i).rep(p).xdata);



            for q = 1:n % for loop to create an array for the fit
               condition(i).rep(p).yfitted(q) = condition(i).rep(p).x(1)...
               *condition(i).rep(p).xdata(q)*...
               (1 + ((2^(1/condition(i).rep(p).x(3))-1) ...
               * condition(i).rep(p).xdata(q)^2)/condition(i).rep(p).x(2))...
               ^-condition(i).rep(p).x(3);
            end
        end       
    end   
end

function [fit, field, data] = corr_time_series(fname)

    if isempty(fname)
        [temp_and_field, data, params] = eprload();
    else
        [temp_and_field, data, params] = eprload(fname);
    end
    field = temp_and_field{1};
    temps = temp_and_field{2};

    % field is the x values which the same for all
    % temps in the temp associated with each spectrum in order
    % each column in data is the y values of the spectrum at the ith
    % temperature step


    %subtract the baseline as defined by the average of the first 100
    %points for all spectra
    for i=1:numel(temps)
        data(:,i) = data(:,i) - mean(data(1:100,i));
    end


    % for the first fit in the middle of the range of the series, allow A,
    % g, and tcorr to vary
    v = floor(numel(temps)/2);

    disp('guessing the starting parameters of the series. This is a full fit, so be patient ...');
    curr_fit = single_component_nitroxide(field, data(:,v), params, []);
    
    rmsd = zeros(numel(temps), 1);
    logtcorr = zeros (numel(temps),3)
    for i=1:numel(temps)
        fit(i) = single_component_nitroxide(field, data(:,i), params, curr_fit);
        curr_fit = fit(i);
        logtcorr(i,1) = curr_fit.argsfit{1,1}.logtcorr(1);
        logtcorr(i,2) = curr_fit.argsfit{1,1}.logtcorr(2);
        rmsd(i) = curr_fit.rmsd;
    end

scatter(temps, logtcorr(:,1), 'Filled');
hold on;

ax = gca;
scatter(ax, temps, logtcorr(:,2), 'Filled');
ax.YLabel.String = 'correlation time (s) (log)';
ax.XLabel.String = 'Temperature (K)';

figure;
ax = gca;
scatter(temps, rmsd, 'Filled');
ax.YLabel.String = 'RMSD';
ax.XLabel.String = 'Temperature (K)';
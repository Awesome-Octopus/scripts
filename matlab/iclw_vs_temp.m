%% desired syntax [linewidths_vs_temp_matrix, graph_handle] = iclw_vs_temp(eprload object_x, eprload object_y, eprload_params)

function [iclw_vs_temp, graph_handle] = iclw_vs_temp(varargin)

    if isempty(varargin)
        [temp_and_field, data, params] = eprload();
    end
    field = temp_and_field{1};
    temps = temp_and_field{2};

    % field is the x values which the same for all
    % temps in the temp associated with each spectrum in order
    % each column in data is the y values of the spectrum at the ith
    % temperature step
    xq = linspace(min(field), max(field), 10000);
    for i=1:numel(temps)

        interpol = griddedInterpolant(field, data(:,i));
        yq = interpol(xq);

        %adjust baseline based on first few points (this assumes that your
        %spectra start with a flat baseline)
        yq = yq - mean(yq(1:100));
        [peaks, interp_field] = findpeaks(yq, xq, 'MinPeakDistance', ...
        10, 'MinPeakProminence', max(yq)/10, 'Annotate',...
        'extents','WidthReference','halfheight');

        %flip the spectrum upside down and find the max
        inverseyq = -1*yq;

        [inversepeaks, inversefield] ...
        = findpeaks(inverseyq, xq, 'MinPeakDistance', ...
        10, 'MinPeakProminence', max(inverseyq)/10, 'Annotate',...
        'extents','WidthReference','halfheight');

        % if the spectrum has 2 or more detectable peaks, the correct peak is the 2nd
        % one
        if numel(peaks) >1 && numel(inversepeaks >1)
            centerpeakx = interp_field(2);
            centertroughx = inversefield(2);
            iclw = 1/abs((centertroughx - centerpeakx));
    
        % otherwise it is the only one
        elseif numel(peaks) < 2 && numel(inversepeaks >1)
             centerpeakx =  field(1);
             centertroughx =  inversefield(1);
             iclw = 1/abs(( centertroughx -  centerpeakx));
    
        else
             warndlg('A spectrum appears to not have an assignable peak, so its width has been set to 1. Go through the code and play with the parameters for the findpeaks function probably.');
              iclw = 1;
        end
        iclws(i) =  iclw;
        t_corr(i) = single_component_nitroxide(field, data(:,i), params)
end

graph_handle = scatter(temps, t_corr, 'Filled');
ax = gca;
ax.YLabel.String = 'logtcor';
ax.XLabel.String = 'Temperature (K)';

set(ax,'Color','None');


    
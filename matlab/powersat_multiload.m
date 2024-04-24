function [condition, doDataMask] = powersat_multiload 
% condition is an indexing variable, so that all experimental conditions 
% can be looped over. 
% condition(1) is nitrogen, condition(2) is air, condition(3) is nickel



filename = uigetfile({'*.dta';'*.spc'},  'Select Nitrogen replicates:','MultiSelect','on');


% if you only select one file it will save filename as a string rather than
% a 1x1 cell array which is what we want
if isa(filename,'cell') == 1
    
    num_nitro_reps = numel(filename);
elseif isa(filename,'char') == 1   
    filename = cellstr(filename);
    num_nitro_reps = 1;
else
    num_nitro_reps = 0;
end
  
if num_nitro_reps ~= 0
    
    %% uncomment this to give a prompt to choose what color to plot in.
    % asks what color to plot nitro in, if none given, black is selected
    %condition(1).color = input('Nitrogen color (default ''k''): ');
    %if isempty(condition(1).color) == 1
    %    condition(1).color = 'k';
    % end
    condition(1).color = 'k';
    
    for i = 1:num_nitro_reps
        [condition(1).rep(i).xdata, condition(1).rep(i).ydata, ...
            condition(1).rep(i).params] = eprload(char(filename(i)));

    end
end
    
filename = uigetfile({'*.dta';'*.spc'},  'Select air replicates:','MultiSelect','on');


% if you only select one file it will save filename as a string rather than
% a 1x1 cell array which is what we want
if isa(filename,'cell') == 1
        
    num_air_reps = numel(filename);        
    
elseif isa(filename,'char') == 1
        filename = cellstr(filename);
        num_air_reps = 1;        
else
        num_air_reps =0;
end

if num_air_reps ~= 0
    
    
    %% uncomment this to give a prompt to choose what color to plot in.
    % asks what color to plot air in, if none given, red is selected
    %condition(2).color = input('Air color (default ''r''): ');
    %if isempty(condition(2).color) == 1
    %    condition(2).color = 'r';
    % end
    condition(2).color = 'r';
    
    for i = 1:num_air_reps
        [condition(2).rep(i).xdata, condition(2).rep(i).ydata, ...
            condition(2).rep(i).params] = eprload(char(filename(i)));
    
    end
end


filename = uigetfile({'*.dta';'*.spc'},  'Select nickel replicates:','MultiSelect','on');


% if you only select one file it will save filename as a string rather than
% a 1x1 cell array which is what we want
if isa(filename,'cell') == 1
    
    num_nickel_reps = numel(filename);
elseif isa(filename,'char') == 1   
    filename = cellstr(filename);
    num_nickel_reps = 1;
else
    num_nickel_reps = 0;
end
       

       
if num_nickel_reps ~= 0
    
    
    %% uncomment this to give a prompt to choose what color to plot in.
    % asks what color to plot nickel in, if none given, blue is selected
    %condition(3).color = input('Nickel color (default ''b''): ');
    %if isempty(condition(3).color) == 1
    %    condition(3).color = 'b';
    % end
    condition(3).color = 'b';
    
    for i = 1:num_nickel_reps
        [condition(3).rep(i).xdata, condition(3).rep(i).ydata, ...
            condition(3).rep(i).params] = eprload(char(filename(i)));

    
    end
end

%% uncomment this to use this script for file in the old winEPR format from the old EMX

% for i = 1:numel(condition)
%    for s = 1:numel(condition(i).rep)
%       condition(i).rep(s).params.YPTS = 17;
%       condition(i).rep(s).params.YMIN = 0.100237;
%       condition(i).rep(s).params.YWID = 158.765409;
%   end
% end

%% set this value to 1 to allow masking of low attenuation points from data analysis
doDataMask = 0;
end
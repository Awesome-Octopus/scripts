% Note: you must have the same start and
% end wavelength for blank as for sample. if one of them extends longer,
% you need to trim it.



[fid , path] = uigetfile('*.dat', 'Select your blank spectrum');
%opens dialog box to select the file

avg_points = input('Perform data smoothing? (enter ''y'' or ''n''): ');
%asks whether 3 point rolling average will be taken

c = input('Convert to mean residue ellipticity?: (''y'' or ''n'') ');
if c == 'y'
    
    do_conversion = 1;
    mrw = str2double(input('Enter mean residue weight (mol. wt./(num. residue - 1)): ', 's'));
    conc = str2double(input('Enter protein concentration in mg/ml: ', 's'));
    
else
    
    do_conversion = 0;
    
end
    
base_adjust = input('Set the baseline to the value at max wavelength? (enter ''y'' or ''n''): ');

fullpath = strcat(path,fid);

fid=fopen(fullpath);
%opens file specified

tline = fgetl(fid);

rawDatalines = cell(0,1);
while ischar(tline)
    rawDatalines{end+1,1} = tline;
    tline = fgetl(fid);
end

%breaks up text read in into 1 line per cell

fclose(fid);

max_wavelength = str2double(cell2mat(regexp(rawDatalines{9},'[0-9]*\.?[0-9]*','match')));
%the 9th line of text contains the metadata for max wavelength.
%this uses a regex to return substrings that are one or more digits
% preceded by ': ', which should be the number in the line

min_wavelength = str2double(cell2mat(regexp(rawDatalines{10},'[0-9]*\.?[0-9]*','match')));
%the 10th line of text contains the metadata for min wavelength.
        
wavelength_step = str2double(cell2mat(regexp(rawDatalines{11},'[0-9]*\.?[0-9]*','match')));
%the 11th line of text contains the metadata for how each measurement is
%spaced apart

num_scans = 1;

rawDatalines = rawDatalines(19:end,1);
%removed first 18 lines which are just metadata, from which we've extracted
%what we need already

 for i = 1:numel(rawDatalines)
%     go through every line
%     
    

    if regexp(rawDatalines{i}, '^\$MDCNAME:Scan') == 1
       
        %the string $MDCNAME:Scan_ is at the start of each of the
        %concatenated scans. count the number of scans by augenting
        %everytime regex matches this string, since it never explicitly
        %lists a number of scans in the header. Notably, this doesnt coun
        %the last averaged data entry, if that is part of the file as put
        %in by during data processing on the instrument
        
            num_scans = num_scans +1;
            
    end

    
    if regexp(rawDatalines{i}, '^(\$|\s)') == 1
    % returns 1 for each line that begins
    % with '$' or a space as these lines
    % are metadata or data labels and need
    % to be removed   
        
        rawDatalines{i} = [];
    end    
    %empty the cell if the regexp gives 1
 end
 
rawDatalines = rawDatalines(~cellfun(@isempty,rawDatalines));
%removed all the emtpy cells 

split_data = regexp(rawDatalines, '[\s]+', 'split');
%breaks cells up into seperate columns whenever it encounters
%any amount of whitespace

trimmedData = cell2table(split_data);
%convert to a table

trimmedData = sortrows(trimmedData);
% this sorts the x y pairs by wavelength ascending so that the first n rows are all
% measurements of the min wavelength where n is the number of scans


fullpath = strcat(fullpath,'_trimmed.txt');

writetable(trimmedData, fullpath);
%the only way to get a proper table with doubles rather than an array of
%cells is to write it to a text file and then reopen the file. if you dont
%convert to a proper table of doubles everything downstream doesnt work


trimmedData = readmatrix(fullpath);



x_range = max_wavelength - min_wavelength;


wavelength = linspace (min_wavelength,max_wavelength,(x_range/wavelength_step+1));
% this assumes data reads every 1 nm not necessarily true, change later

for i = 1:numel(wavelength)
% go through each wavelength
  
    for p=1:num_scans
    % in the sorted data array, each block of x rows are measurements from 
    % different scans at the same wavelength
        
        if (i == 1)
           
            CD_signal(1,p) = p;
            CD_signal(1,p) = trimmedData(p,2);
        % this silly if statement has to be here because it doesn't work if i = 1
        % and matlab doesn't allow arrays to start at 0, like literally
        % every other programming language, which would fix this.
        
    
        else
            
            CD_signal(i,p) = (i-1)*num_scans + p;
            CD_signal(i,p) = trimmedData((i-1)*num_scans + p,2);
            %make a new array of vectors containing the readings at each
            %wavelength
            
        
        end
        
       %mean_CD(i) = mean(CD_signal(i,:)); 
       %stdev_CD(i) = std(CD_signal(i,:)); 
    end  

%     stdev_CD = std(CD_signal,2);
end


CD_signal = zeros (numel(wavelength),num_scans);
blank_mean_CD = zeros (numel(wavelength),1);
%stdev_CD = zeros (x_range,1);
%pre-initialling variables for speed


% this should be rewritten using a modulo operator so that every X rows it
% moves to a new wavelength index and restarts summation and averaging

for i = 1:numel(wavelength)
%go through each wavelength
  
    for p=1:num_scans
    % in the sorted data array, each block of x rows are measurements from 
    %different scans at the same wavelength
        
        if (i == 1)
            
            CD_signal(1,p) = p;
            CD_signal(1,p) = trimmedData(p,2);
        % this silly if statement has to be here because it doesn't work if i = 1
        % and matlab doesn't allow arrays to start at 0, like literally
        % every other programming language, which would fix this.
        
    
        else
            
            CD_signal(i,p) = (i-1)*num_scans + p;
            CD_signal(i,p) = trimmedData((i-1)*num_scans + p,2);
            %make a new array of vectors containing in each row the readings at each
            %from each scan. each row is a wavelength, each column is a
            %reading from a scan at that wavelength
            
        
        end
        
       blank_mean_CD(i) = mean(CD_signal(i,:)); 
       % blank_stdev_CD(i) = std(CD_signal(i,:)); 
    end  

end

if avg_points == 'y'
    
    mean_CD = datasmooth(blank_mean_CD, 1);
    % this performs a rolling average over the mean data 
end

%% I'm doing this the inelegant way and copy pasting code rather than define functions and call them
% this overwites all the previous variables except for blank_mean_CD and blank_stdev_CD
% and min and max wavelength

[fid , path] = uigetfile('*.dat', 'Select your sample spectrum');
%opens dialog box to select the file

fullpath = strcat(path,fid);


fid=fopen(fullpath);
%opens file specified

tline = fgetl(fid);

rawDatalines = cell(0,1);
while ischar(tline)
    rawDatalines{end+1,1} = tline;
    tline = fgetl(fid);
end

%breaks up text read in into 1 line per cell

fclose(fid);

rawDatalines = rawDatalines(19:end,1);
%removed first 18 lines which are just metadata

num_scans = 1;

 for i = 1:numel(rawDatalines)
%     go through every line
%     
    

    if regexp(rawDatalines{i}, '^\$MDCNAME') == 1
       
        %the string $MDCNAME:Scan_ is at the start of each of the
        %concatenated scans. count the number of scans by augenting
        %everytime regex matches this string, since it never explicitly
        %lists a number of scans in the header
        
            num_scans = num_scans +1;
            
    end
 end
 for i = 1:numel(rawDatalines)
%     go through every line
%     
    
    rex = regexp(rawDatalines{i}, '^(\$|\s)');
    % returns 1 for each line that begins
    % with '$' or a space as these lines
    % are metadata or data labels and need
    % to be removed
           
    if rex == 1
       rawDatalines{i} = [];
    end    
    %empty the cell if the regexp gives 1
 end
 
rawDatalines = rawDatalines(~cellfun(@isempty,rawDatalines));
%removed all the emtpy cells 

split_data = regexp(rawDatalines, '[\s]+', 'split');
%breaks cells up into seperate columns whenever it encounters
%any amount of whitespace

trimmedData = cell2table(split_data);
%convert to a table

trimmedData = sortrows(trimmedData);

fullpath = strcat(fullpath,'_result.txt');

writetable(trimmedData, fullpath, 'Delimiter', 'tab');
%the only way to get a proper table with doubles rather than an array of
%cells is to write it to a text file and then reopen the file. if you dont
%convert to a proper table of doubles everything downstream doesnt work


trimmedData = readmatrix(fullpath);


for i = 1:x_range
% go through each wavelength
  
    for p=1:num_scans
    % in the sorted data array, each block of x rows are measurements from 
    % different scans at the same wavelength
        
        if (i == 1)
           
            CD_signal(1,p) = p;
            CD_signal(1,p) = trimmedData(p,2);
        % this silly if statement has to be here because it doesn't work if i = 1
        % and matlab doesn't allow arrays to start at 0, like literally
        % every other programming language, which would fix this.
        
    
        else
            
            CD_signal(i,p) = (i-1)*num_scans + p;
             CD_signal(i,p) = trimmedData((i-1)*num_scans + p,2);
            %make a new array of vectors containing the readings at each
            %wavelength
            
        
        end
        
       %mean_CD(i) = mean(CD_signal(i,:)); 
       %stdev_CD(i) = std(CD_signal(i,:)); 
    end  

%     stdev_CD = std(CD_signal,2);
end


CD_signal = zeros (numel(wavelength),num_scans);
mean_CD = zeros (numel(wavelength),1);
stdev_CD = zeros (numel(wavelength),1);
%pre-initialling variables for speed

for i = 1:numel(wavelength)
%go through each wavelength
  
    for p=1:num_scans
    % in the sorted data array, each block of x rows are measurements from 
    %different scans at the same wavelength
        
        if (i == 1)
            
            CD_signal(1,p) = p;
            CD_signal(1,p) = trimmedData(p,2);
        % this silly if statement has to be here because it doesn't work if i = 1
        % and matlab doesn't allow arrays to start at 0, like literally
        % every other programming language, which would fix this.
        
    
        else
            
            CD_signal(i,p) = (i-1)*num_scans + p;
            CD_signal(i,p) = trimmedData((i-1)*num_scans + p,2);
            %make a new array of vectors containing in each row the readings at each
            %from each scan. each row is a wavelength, each column is a
            %reading from a scan at that wavelength
            
        
        end
        
       mean_CD(i) = mean(CD_signal(i,:)); 
       %stdev_CD(i) = std(CD_signal(i,:)); 
    end  

end

if avg_points == 'y'
    
    mean_CD = datasmooth(mean_CD, 1);
    
end

CD_result(:,1) = wavelength;
CD_result(:,2) = (mean_CD - blank_mean_CD);

if base_adjust == 'y'
    
    CD_result(:,2) = CD_result(:,2) - CD_result(numel(wavelength),2);
end

if do_conversion == 1
    
    CD_result(:,2) = CD_result(:,2)*mrw/conc;
end    

plot (CD_result(:,1), CD_result(:,2));

if do_conversion == 1
    
    ylabel('mean residue ellipticity (mdeg*cm^2*dmol^-^1)')
else
   
    ylabel('ellipticity (mdeg)');
    
end
xlabel('wavelength (nm)');

writematrix(CD_result, fullpath, 'Delimiter', 'tab');
%writes a csv for the buffer subtracted spectrum.

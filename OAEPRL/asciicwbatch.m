%% CW EPR Batch ACSII Conversion Tool
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.14
% Last Modified 09/29/2015

clear all

%% Creates a structure named files with a list of all the files in the
% current directory

path = uigetdir()

%% A for loop that cycles through all of the files in the directory, loads
% them in using EasySpin, creates a filename with a .txt extension from the
% original filename and saves an ASCII file

pathextent = sprintf('%s/*.spc',path);

files = dir(pathextent);

L = length(files);

for i = 1:L % runs a for loop through the files
    [x,y] = eprload(files(i).name);
    fname = files(i).name(1:end-4) % removes the *spc from the filename
    textname = sprintf('%s.txt',fname); % adds .txt
    data = [x',y'];
    save(textname,'data','-ascii');
end


pathextent = sprintf('%s/*.DTA',path);

files = dir(pathextent);

L = length(files);

for i = 1:L % runs a for loop through the files
    [x,y] = eprload(files(i).name);
    fname = files(i).name(1:end-4) % removes the *spc from the filename
    textname = sprintf('%s.txt',fname); % adds .txt
    data = [x,y];
    save(textname,'data','-ascii');
end
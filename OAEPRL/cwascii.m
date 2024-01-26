%% CW ASCII Directory
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

clear all; % clears all of the variables and structures

clear all

files = dir('*.spc');
[p q] = size(files);

for i=1:p
    
    [x y params fname] = eprload(files(i).name);
    fname = files(i).name(1:end-4)
    fname = sprintf('%s.txt',fname)
    data = [x',y']
    save(fname,'data','-ascii')
end
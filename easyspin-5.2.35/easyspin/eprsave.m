% eprsave  Save data in Bruker BES3T format
%
%   eprsave(FileName,x,data)
%   eprsave(FileName,x,data,TitleString)
%   eprsave(FileName,x,data,TitleString,mwFreq)
%
%   Saves the dataset in x and data in the Bruker BES3T
%   format in a .DTA and a .DSC file with the file
%   name given in FileName. x is the x axis data, and data
%   is the intensity data (real or complex). TitleString is
%   the name of the dataset that will be displayed in the
%   Bruker software. mwFreq is the microwave frequency, in GHz.
%
%   Two-dimensional data can be saved by giving a matrix in data,
%   and by supplying both axes in x as a cell array.
%
%   Examples:
%     eprsave(myFilename,B,spc);       % save 1D data spc, x axis = B 
%     eprsave(myFilename,{t1,t2},V);   % save 2D data matrix V,
%                                      %   x axis = t1, y axis = t2

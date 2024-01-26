% levelsplot    plot energy level (Zeeman) diagram 
%
%  levelsplot(Sys,Ori,B)
%  levelsplot(Sys,Ori,B,mwFreq)
%  levelsplot(Sys,Ori,B,mwFreq,Par)
%
%    Sys        spin system structure
%    Ori        orientation of magnetic field in molecular frame
%               - string: 'x','y','z','xy','xz','yz', or 'xyz'
%               - 2-element vector [phi theta]
%               - 3-element vector [phi theta chi]
%    B          field range, in mT; either Bmax, [Bmin Bmax], or a full vector
%    mwFreq     spectrometer frequency, in GHz
%    Par        other parameters
%      Units           energy units for plotting, 'GHz' or 'cm^-1' or 'eV'
%      nPoints         Number of points
%      ColorThreshold  Coloring threshold. All transitions with relative
%                      intensity below this will be gray. Example: 0.05
%      PlotThreshold   All transitions below with relative intensity
%                      below this value will not be plotted. Example: 0.005
%
%  If mwFreq is given, resonances are drawn. Red lines indicate allowed
%  transitions, gray lines forbidden ones. If the lines are terminated with
%  dots, the relative transition intensity is larger than 1%.
%
%  Example:
%    Sys = struct('S',7/2,'g',2,'D',5e3);
%    levelsplot(Sys,'xy',6e3,95);

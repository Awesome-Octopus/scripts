% levelsplot    Plot energy level (Zeeman) diagram 
%
%  levelsplot(Sys,Ori,B)
%  levelsplot(Sys,Ori,B,mwFreq)
%  levelsplot(Sys,Ori,B,mwFreq,Opt)
%
%    Sys        spin system structure
%    Ori        (a) orientation of magnetic field vector in molecular frame
%               - string: 'x','y','z','xy','xz','yz', or 'xyz'
%               - 2-element vector [phi theta] (radians)
%               (b) orientation of lab frame in molecular frame
%               - 3-element vector [phi theta chi] (radians)
%    B          field range, in mT; either Bmax, [Bmin Bmax], or a full vector
%    mwFreq     spectrometer frequency, in GHz
%    Opt        options
%      Units           energy units for plotting, 'GHz' or 'cm^-1' or 'eV'
%      nPoints         number of points
%      PlotThreshold   all transitions with relative intensity below
%                      this value will not be plotted. Example: 0.005
%      SlopeColor      true/false (default false). Color energy level lines
%                      by their slope, (corresponding to their mS expectation
%                      value).
%      StickSpectrum   true/false (default false). Plot stick spectrum underneath
%                      Zeeman diagram. Default is false.
%
%  If mwFreq is given, resonances are drawn. Red lines indicate allowed
%  transitions, gray lines forbidden ones. Hovering with the cursor over
%  the lines displays intensity information.
%
%  Example:
%    Sys = struct('S',7/2,'g',2,'D',5000);
%    levelsplot(Sys,'xy',[0 6000],95);

% stackplot    Stacked plot of data
%
%    stackplot(x,y)
%    stackplot(x,y,scale)
%    stackplot(x,y,scale,step)
%    stackplot(x,y,scale,step,sliceLabels)
%    stackplot(x,y,scale,step,sliceLabels,colors)
%
%  Plots the columns or rows of y stacked on top of each other. The length
%  of x determines whether columns or rows are plotted. 
%  Slices are rescaled as specified in scale ('maxabs', 'int', 'dint', or 
%  'none', see rescaledata.m for details).
%
% Input: x     - vector of x-axis values
%        y     - matrix of data to plot
%        scale - string defining the type of scaling, options are 'maxabs',
%                'int' (normalized integral), 'dint' (normalized double integral), 
%                or 'none' (default = 'maxabs')
%                an additional scaling factor can be defined by providing a
%                cell as input {scale scalefact}, e.g. {'maxabs',5} (useful if 
%                the shifts in the stacked plot correspond to a specific axis)
%        step  - relative step size for the stacked plot (as a fraction of
%                the maximum amplitude of the plotted data) or list of positions
%                for the different slices on the y-axis
%                (default = 1.5)
%        sliceLabels - list of values or cell array for labeling of the
%                      slices on the y-axis
%        colors      - single color, list of colors or string with colormap name
%                      for plotting of the different slices
%

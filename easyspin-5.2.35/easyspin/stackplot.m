% stackplot    Stacked plot of data
%
%    stackplot(x,y)
%    stackplot(x,y,scale)
%    stackplot(x,y,scale,step)
%
%    Plots the columns or rows of y stacked on
%    top of each other. The length of x determines
%    whether columns or rows are plotted. Slices
%    are rescaled to max-min=scale and plotted
%    at distance of step.
%
%    If step is missing, step = 1 is assumed.
%    If scale is missing, scale = 1 is assumed.
%    If scale is set to 0, no rescaling is done.

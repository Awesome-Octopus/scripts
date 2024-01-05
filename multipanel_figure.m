% this will make a n by m multipanel figure out of independent plots each
% of which is saved as .fig in the current directory. select the figures 
% from the ui in the order you want them to appear left to right row-by row

%% input

% %open a dialog box that gives instructions for the first prompt
% instructions = uifigure;
% uialert(instructions, "Select the figures in the order you want them to appear, left to right, top to bottom",...
% "Instructions", 'Icon', 'info');

%prompt the user for how many rows and how many columns
response = inputdlg({'Number of rows of panels:', 'Number of columns of panels:'}, ...
'Enter Dimensions', [1 90]);

r = str2double(response{1});
c = str2double(response{2});
filenames = cell(r,c);
axesArray = cell(r,c);

for i=1:r
    
   for p=1:c
               
       
       % asks user to select a file for each panel in the order they are to appear from a file selction box
       filenames{i,p} = uigetfile("*.fig", sprintf(...
       "select panel for row #%d, column #%d", i, p),"MultiSelect", "off");   
   

   
       f = openfig(filenames{i,p}, 'invisible');
       a = gca;
%      a.YLabel.String = 'Norm Ampl';
       a.YLabel.FontSize = 5;
%      a.YLimMode = 'Manual';
%      a.YLim = [0, 1.2];
       delete(findobj(a,'Type', 'Text'));
       axesArray{i,p} = a;
    end    
end
 
 figure;
 tcl = tiledlayout(r,c);
 tileNum = 0;

for i=1:6
    
   for p=1:3    
       tileNum = tileNum +1;
       axesArray{i,p}.FontSize = 4;
       gca.Children.MarkerSize = 1;
       axesArray{i,p}.Parent = tcl;
       axesArray{i,p}.Layout.Tile = tileNum;             
    end    
end
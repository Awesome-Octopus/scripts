%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  W  A  R  N  I  N  G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% currently this program will not work on MATLAB r2023a, it works best on
% r2020a

%% a GUI to visualize distance distributions for proteins from MD trajectories
%  this should be used in conjuction with the scripts "write_xyz.tcl" and
%  "write_xyz.slurm", which runs a script en batch on VMD that prints out the
%  backbone positions and side chain positions for every amino acid in 
%  a specified range of frames.
%  You may need to read through and modify those scripts to have the
%  necessary input to run this script.

%% still needed:
% 1. find a way to calculate the optimum binning number for all label pairs,
% and then force that binning number for all the histograms so that they
% have the same area under the curve
% 2. implement TOAC label distributions based on a static position
% relative to the backbone
% 3. implement a way to remove label pairs from GUI
% 4. implement MTSL label distribution based on simple extension of the
% vector between BB and SC1 or between alpha and beta carbon
% 5. implement MTSL label distribution based on a cone model of rotamers, 
% possibly with arbitrary scaling of the flexibility of the label at given
% sites.
% 6. integrate the running of the of the VMD script with this to be more
% user friendly, possibly with GUI tools to allow for the selection of
% different sets of sites for writing the d.dat file to save time.


%% import
clear all;
w = readtable(uigetfile("*.dat"),"FileType","text","Delimiter","\t",...
    "ReadVariableNames",true,"NumHeaderLines",1);

% depending on the default import functions of what version of matlab we
% are using we may get the column "frame_" as the first and a junk column
% as the last. if we do we need to trim these off
if strcmp(w.Properties.VariableNames{1}, 'frame_') == 1
    w = w(:,2:end);
end

% the last data column should end in '_Z', if not it is likely an untrimmed
% column of newline characters that was imported as data and should be
% erased
while strcmp(w.Properties.VariableNames{end}(end), 'Z') == 0
    w = w(:,1:end-1);
end


sitenames = w.Properties.VariableNames;
data = table2array(w);

%%
clear w;
[numframes, m] = size(data);

% since each site has an x y and z coordinate
numsites = m/3;

% for debugging, if there is not x y and z for each site 
% numsites != an integer, let user know there is a formatting problem
if mod(numsites, 1) ~= 0
    fprintf("number of coordinate columns is not evenly divisible by 3. Check input formatting or change readtable options on line 35")
end    

% if any vectors in a frame are NaN, throw the whole frame out
has_nans = false;
for i=numframes
    if sum(isnan(data(i,:))) ~= 0
        data(i,:) = [];
        numframes = numframes - 1;
        has_nans = true;
    end
end    

position = cell(numframes,numsites);

% read in each row as sets of three columns at a time which define
% a spatial vector
for i=1:numframes

    q = 1;
    for j=1:3:m
        position(i,q) = {[data(i,j), data(i,j+1), data(i,j+2)]};
        q = q + 1;
    end

end
clear q;

clear data;


% one name for each x y z column
sitenames = sitenames(1:3:end);

% the base figure which we will be working off of
fig = uifigure("Position", [50, 50, 1100, 500],"UserData", ...
    struct( "guiobjs", []), "Visible", "on", "Name", "MD2DEER");


% assign identifiers to the positions based on column headers
% this will contain all position and sequence data
% organized as chain -> resnum -> BB or side chain -> coordinates
coordinate = struct('name',strings,'type', strings,'seq', strings, ...
    'chain', strings, 'xyz', zeros(1,3));
for i=1:numel(sitenames)
    str = regexp(sitenames{i}, "_", 'split');

    coordinate(i).name = str{1};
    coordinate(i).type = str{2};
    coordinate(i).seq = str{3};
    coordinate(i).chain = str{4};
    coordinate(i).xyz = position(:,i);

end
clear str;
% make a list that the user can select from for the chain letter and 
% sequence number based on tags pulled from column headers 
seqlist(1) = {coordinate(1).seq};
chainlist(1) = {coordinate(1).chain};
for i = 2:numel(coordinate)
    if strcmpi(coordinate(i).seq, coordinate(i-1).seq) == false
        seqlist(end+1) = {coordinate(i).seq};
    end    
    if strcmpi(coordinate(i).chain, coordinate(i-1).chain) == false
        chainlist(end+1) = {coordinate(i).chain};
    end    
end    

chainlist = unique(chainlist, "stable");
seqlist = unique(seqlist, "stable");

%% set up ui


% ask the user which sites he wants

addsitebtn = uibutton(fig, 'push', 'Position', [100, 450, 125, 22], ...
    'Text', 'Add Label Pair', 'ButtonPushedFcn', ...
    @(addsitebtn, event) addsite(addsitebtn, fig, seqlist, chainlist), ...
    'Tooltip', 'Add selection for additional label pair');

calculatebtn = uibutton(fig, 'push', 'Position', [300 450 125 22], ...
    'Text', 'Get Distances', 'ButtonPushedFcn', ...
    @(calculatebtn, event) update(fig, coordinate,...
    numframes));

% pick a label model type
labellist = {'MTSL' 'TOAC' 'Backbone'};
labeltype = uidropdown (fig, "Items", labellist, ...
        "ItemsData", labellist, "Position", ...
        [120 20 90 22]);
lbl1 = uilabel(fig, "Position", [130 40 130 22]);
lbl1.Text = "Label Type";

ShowComponentSwitch = uiswitch(fig, 'slider', 'position', [370 20 100 22], ...
    'ItemsData', [0 1]);
lbl2 = uilabel(fig, "Position", [350 40 130 22]);
lbl2.Text = "Show Components";


% display plot
ax = uiaxes('Parent',fig,'Position',[500, 40, 560, 420]);
ax.XLabel.String = "Distance ({\AA})";
ax.XLabel.Interpreter = "latex";
ax.YLabel.String = "Frequency";
ax.XLimMode = "manual";
ax.XLim = [0 80];
%% OAEPRL Processing Shortcuts Generator 1/11/2011 1.13
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

clear all; % clears all of the variables and structures

iconfile = which('oaeprl.gif'); % Finds the icon in the Matlab paths
isEditable = 'true'; % Allows the shortcuts to be writeable

name = 'ASCII Conversion'; % The shortcut name
evalstr = 'textconvert'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut


name = 'ASCII Batch Conversion'; % The shortcut name
evalstr = 'asciiwbatch'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut


name = 'ASCII Plotting'; % The shortcut name
evalstr = 'asciidata'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut

name = 'Concentration'; % The shortcut name
evalstr = 'concentration'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut

name = 'Comparison'; % The shortcut name
evalstr = 'comparison'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut

name = 'Power Saturation'; % The shortcut name
evalstr = 'powersat'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut

name = 'Power Saturation Batch'; % The shortcut name
evalstr = 'powersatbatch'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut

name = 'T1 Fit'; % The shortcut name
evalstr = 't1fit'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut

name = 'T2 Fit'; % The shortcut name
evalstr = 't2fit'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut

name = 'T2 Fit-DEER'; % The shortcut name
evalstr = 't2fit_DEER'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut

name = '2P ESEEM'; % The shortcut name
evalstr = 'ESEEM2pproc'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut

name = '3P ESEEM'; % The shortcut name
evalstr = 'ESEEM3pproc'; % Shortcut target

awtinvoke(com.mathworks.mlwidgets.shortcuts.ShortcutUtils,'addShortcutToBottom',name,evalstr,iconfile,'Toolbar Shortcuts',isEditable); % Creates the shortcut


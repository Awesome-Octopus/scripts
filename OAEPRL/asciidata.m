%% ASCII Data Plotter
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

function varargout = asciidata(varargin)
% ASCIIDATA M-file for asciidata.fig
%      ASCIIDATA, by itself, creates a new ASCIIDATA or raises the existing
%      singleton*.
%
%      H = ASCIIDATA returns the handle to a new ASCIIDATA or the handle to
%      the existing singleton*.
%
%      ASCIIDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASCIIDATA.M with the given input arguments.
%
%      ASCIIDATA('Property','Value',...) creates a new ASCIIDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before elexsysdata_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to asciidata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help asciidata

% Last Modified 02/21/2011 GUIDE v2.5 17-Feb-2011 16:31:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @asciidata_OpeningFcn, ...
                   'gui_OutputFcn',  @asciidata_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before asciidata is made visible.
function asciidata_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to asciidata (see VARARGIN)

% Choose default command line output for asciidata
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes asciidata wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = asciidata_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in openA.
function openA_Callback(hObject, eventdata, handles)
% hObject    handle to openA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[colA1 colA2 colA3 colA4 colA5 colA6 colA7 colA8 colA9] = textread(uigetfile('*.txt'),'%f%f%f%f%f%f%f%f%f','headerlines',2);
dataA = [colA1 colA2 colA3 colA4 colA5 colA6 colA7 colA8 colA9];
handles.dataAx= dataA(:,handles.colAx);
handles.dataAy= dataA(:,handles.colAy);

guidata(hObject, handles)

% --- Executes on button press in openB.
function openB_Callback(hObject, eventdata, handles)
% hObject    handle to openB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[colB1 colB2 colB3 colB4 colB5 colB6 colB7 colB8 colB9] = textread(uigetfile('*.txt'),'%f%f%f%f%f%f%f%f%f','headerlines',2);
dataB = [colB1 colB2 colB3 colB4 colB5 colB6 colB7 colB8 colB9];
handles.dataBx= dataB(:,handles.colBx);
handles.dataBy= dataB(:,handles.colBy);

guidata(hObject, handles)

function Ax_Callback(hObject, eventdata, handles)
% hObject    handle to Ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ax as text
%        str2double(get(hObject,'String')) returns contents of Ax as a double

colAx = str2double(get(hObject, 'String'));
handles.colAx = colAx;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Ax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ay_Callback(hObject, eventdata, handles)
% hObject    handle to Ay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ay as text
%        str2double(get(hObject,'String')) returns contents of Ay as a double

colAy = str2double(get(hObject, 'String'));
handles.colAy = colAy;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Ay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function By_Callback(hObject, eventdata, handles)
% hObject    handle to By (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of By as text
%        str2double(get(hObject,'String')) returns contents of By as a double

colBy = str2double(get(hObject, 'String'));
handles.colBy = colBy;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function By_CreateFcn(hObject, eventdata, handles)
% hObject    handle to By (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Bx_Callback(hObject, eventdata, handles)
% hObject    handle to Bx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bx as text
%        str2double(get(hObject,'String')) returns contents of Bx as a double

colBx = str2double(get(hObject, 'String'));
handles.colBx = colBx;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Bx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dataAadd_Callback(hObject, eventdata, handles)
% hObject    handle to dataAadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataAadd as text
%        str2double(get(hObject,'String')) returns contents of dataAadd as a double

addA = str2double(get(hObject, 'String'));
handles.addA = addA;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function dataAadd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataAadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dataAsub_Callback(hObject, eventdata, handles)
% hObject    handle to data1sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data1sub as text
%        str2double(get(hObject,'String')) returns contents of data1sub as a double

subA = str2double(get(hObject, 'String'));
handles.subA = subA;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function dataAsub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data1sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dataAmult_Callback(hObject, eventdata, handles)
% hObject    handle to dataAmult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataAmult as text
%        str2double(get(hObject,'String')) returns contents of dataAmult as a double

multA = str2double(get(hObject, 'String'));
handles.multA = multA;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function dataAmult_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataAmult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dataAdiv_Callback(hObject, eventdata, handles)
% hObject    handle to dataAdiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataAdiv as text
%        str2double(get(hObject,'String')) returns contents of dataAdiv as a double

divA = str2double(get(hObject, 'String'));
handles.divA = divA;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function dataAdiv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataAdiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dataBadd_Callback(hObject, eventdata, handles)
% hObject    handle to dataBadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataBadd as text
%        str2double(get(hObject,'String')) returns contents of dataBadd as a double

addB = str2double(get(hObject, 'String'));
handles.addB = addB;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function dataBadd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataBadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dataBsub_Callback(hObject, eventdata, handles)
% hObject    handle to dataAsub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataAsub as text
%        str2double(get(hObject,'String')) returns contents of dataAsub as a double

subB = str2double(get(hObject, 'String'));
handles.subB = subB;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.

function dataBsub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataAsub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dataBmult_Callback(hObject, eventdata, handles)
% hObject    handle to dataBmult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataBmult as text
%        str2double(get(hObject,'String')) returns contents of dataBmult as a double

multB = str2double(get(hObject, 'String'));
handles.multB = multB;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function dataBmult_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataBmult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dataBdiv_Callback(hObject, eventdata, handles)
% hObject    handle to dataBdiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataBdiv as text
%        str2double(get(hObject,'String')) returns contents of dataBdiv as a double

divB = str2double(get(hObject, 'String'));
handles.divB = divB;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function dataBdiv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataBdiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addA.
function addA_Callback(hObject, eventdata, handles)
% hObject    handle to addA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataAy = handles.dataAy + handles.addA;
handles.dataAy = dataAy;
plot(handles.dataAx,handles.dataAy,'b',handles.dataBx,handles.dataBy,'r');
guidata(hObject,handles)

% --- Executes on button press in subA.
function subA_Callback(hObject, eventdata, handles)
% hObject    handle to subA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataAy = handles.dataAy-handles.subA;
handles.dataAy = dataAy;
plot(handles.dataAx,handles.dataAy,'b',handles.dataBx,handles.dataBy,'r');
guidata(hObject,handles)

% --- Executes on button press in multA.
function multA_Callback(hObject, eventdata, handles)
% hObject    handle to multA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataAy = handles.dataAy.*handles.multA;
handles.dataAy = dataAy;
plot(handles.dataAx,handles.dataAy,'b',handles.dataBx,handles.dataBy,'r');
guidata(hObject,handles)

% --- Executes on button press in divA.
function divA_Callback(hObject, eventdata, handles)
% hObject    handle to divA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataAy = handles.dataAy./handles.divA;
handles.dataAy = dataAy;
plot(handles.dataAx,handles.dataAy,'b',handles.dataBx,handles.dataBy,'r');
guidata(hObject,handles)

% --- Executes on button press in addB.
function addB_Callback(hObject, eventdata, handles)
% hObject    handle to addB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataBy = handles.dataBy+handles.addB;
handles.dataBy = dataBy;
plot(handles.dataAx,handles.dataAy,'b',handles.dataBx,handles.dataBy,'r');
guidata(hObject,handles)

% --- Executes on button press in subB.
function subB_Callback(hObject, eventdata, handles)
% hObject    handle to subB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataBy = handles.dataBy-handles.subB;
handles.dataBy = dataBy;
plot(handles.dataAx,handles.dataAy,'b',handles.dataBx,handles.dataBy,'r');
guidata(hObject,handles)

% --- Executes on button press in multB.
function multB_Callback(hObject, eventdata, handles)
% hObject    handle to multB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataBy = handles.dataBy.*handles.multB;
handles.dataBy = dataBy;
plot(handles.dataAx,handles.dataAy,'b',handles.dataBx,handles.dataBy,'r');
guidata(hObject,handles)

% --- Executes on button press in divB.
function divB_Callback(hObject, eventdata, handles)
% hObject    handle to divB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataBy = handles.dataBy./handles.divB;
handles.dataBy = dataBy;
plot(handles.dataAx,handles.dataAy,'b',handles.dataBx,handles.dataBy,'r');
guidata(hObject,handles)

function filenameB_Callback(hObject, eventdata, handles)
% hObject    handle to filename2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename2 as text
%        str2double(get(hObject,'String')) returns contents of filename2 as a double

filenameB = get(hObject, 'String');
handles.filenameB = filenameB;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function filenameB_CreateFcn(hObject, eventdata, handles)

% hObject    handle to filename2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot(handles.dataAx,handles.dataAy,'r',handles.dataBx,handles.dataBy,'b');

% --- Executes on button press in writeascii.
function writeascii_Callback(hObject, eventdata, handles)
% hObject    handle to writeascii (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

combdataA = [handles.dataAx handles.dataAy];
combdataB = [handles.dataBx handles.dataBy];
dlmwrite(handles.filenameA,combdataA,'\t');
dlmwrite(handles.filenameB,combdataB,'\t');




function filenameA_Callback(hObject, eventdata, handles)
% hObject    handle to filenameA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filenameA as text
%        str2double(get(hObject,'String')) returns contents of filenameA as a double

filenameA = get(hObject, 'String');
handles.filenameA = filenameA;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function filenameA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filenameA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1

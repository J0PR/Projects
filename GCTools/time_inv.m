function varargout = time_inv(varargin)
% TIME_INV MATLAB code for time_inv.fig
%      TIME_INV, by itself, creates a new TIME_INV or raises the existing
%      singleton*.
%
%      H = TIME_INV returns the handle to a new TIME_INV or the handle to
%      the existing singleton*.
%
%      TIME_INV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIME_INV.M with the given input arguments.
%
%      TIME_INV('Property','Value',...) creates a new TIME_INV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before time_inv_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to time_inv_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help time_inv

% Last Modified by GUIDE v2.5 18-Jun-2015 22:52:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @time_inv_OpeningFcn, ...
                   'gui_OutputFcn',  @time_inv_OutputFcn, ...
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


% --- Executes just before time_inv is made visible.
function time_inv_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to time_inv (see VARARGIN)

% Choose default command line output for time_inv
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes time_inv wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = time_inv_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in cbx_model.
function cbx_model_Callback(hObject, eventdata, handles)
% hObject    handle to cbx_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cbx_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cbx_model


% --- Executes during object creation, after setting all properties.
function cbx_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cbx_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function varargout = sid_temp(varargin)
% SID_TEMP Application M-file for sid_temp.fig
%   SID_TEMP, by itself, creates a new SID_TEMP or raises the existing
%   singleton*.
%
%   H = SID_TEMP returns the handle to a new SID_TEMP or the handle to
%   the existing singleton*.
%
%   SID_TEMP('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in SID_TEMP.M with the given input arguments.
%
%   SID_TEMP('Property','Value',...) creates a new SID_TEMP or raises the
%   existing singleton*.  Starting from the left, property value pairs are
%   applied to the GUI before lbox2_OpeningFunction gets called.  An
%   unrecognized property name or invalid value makes property application
%   stop.  All inputs are passed to sid_temp_OpeningFcn via varargin.
%
%   *See GUI Options - GUI allows only one instance to run (singleton).
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2000-2006 The MathWorks, Inc.

% Edit the above text to modify the response to help sid_temp

% Last Modified by GUIDE v2.5 14-Mar-2012 15:00:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',          mfilename, ...
                   'gui_Singleton',     gui_Singleton, ...
                   'gui_OpeningFcn',    @sid_temp_OpeningFcn, ...
                   'gui_OutputFcn',     @sid_temp_OutputFcn, ...
                   'gui_LayoutFcn',     [], ...
                   'gui_Callback',      []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before sid_temp is made visible.
function sid_temp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sid_temp (see VARARGIN)

% Choose default command line output for sid_temp

addpath ./HspiceToolbox;
handles.output = hObject;

% Update handles structure
handles.identity = 0;
handles.current_choice = 'N4SID';
guidata(hObject, handles);

if nargin == 3,
    initial_dir = pwd;
elseif nargin > 4
    if strcmpi(varargin{1},'dir')
        if exist(varargin{2},'dir')
            initial_dir = varargin{2};
        else
            errordlg('Input argument must be a valid directory','Input Argument Error!')
            return
        end
    else
        errordlg('Unrecognized input argument','Input Argument Error!');
        return;
    end
end
% Populate the listbox
load_listbox(initial_dir,handles)
% Return figure handle as first output argument
    
% UIWAIT makes sid_temp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sid_temp_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ------------------------------------------------------------
% Callback for list box - open .fig with guide, otherwise use open
% ------------------------------------------------------------
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

get(handles.figure1,'SelectionType');
if strcmp(get(handles.figure1,'SelectionType'),'open')
    index_selected = get(handles.listbox1,'Value');
    file_list = get(handles.listbox1,'String');
    filename = file_list{index_selected};
    if  handles.is_dir(handles.sorted_index(index_selected))
        cd (filename)
        load_listbox(pwd,handles)
    else
        [path,name,ext,ver] = fileparts(filename);
        switch ext
            case '.fig'
                guide (filename)
            case '.mat'
                load(filename);
                handles.y_temp_vec = y_temp_vec;
                handles.u_p_vec = u_p_vec;
                handles.time = t;
                guidata(hObject, handles);
            otherwise
                try
                    open(filename)
                catch
                    errordlg(lasterr,'File Type Error','modal')
                end
        end
    end
end
% ------------------------------------------------------------
% Read the current directory and sort the names
% ------------------------------------------------------------
function load_listbox(dir_path,handles)
cd (dir_path)
dir_struct = dir(dir_path);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.file_names = sorted_names;


handles.is_dir = [dir_struct.isdir];
handles.sorted_index = sorted_index;
guidata(handles.figure1,handles)
set(handles.listbox1,'String',handles.file_names,...
	'Value',1)
%set(handles.text8,'String',pwd)


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Add the current directory to the path, as the pwd might change thru' the
% gui. Remove the directory from the path when gui is closed 
% (See figure1_DeleteFcn)
setappdata(hObject, 'StartPath', pwd);
addpath(pwd);
addpath ./HspiceToolbox;



% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Remove the directory added to the path in the figure1_CreateFcn.
if isappdata(hObject, 'StartPath')
    rmpath(getappdata(hObject, 'StartPath'));
end



% --- Executes on button press in SID_button.
function SID_button_Callback(hObject, eventdata, handles)
% hObject    handle to SID_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
%display('SID');
app_start = 1;

try
    order = handles.order;
    app_end = handles.app_end;
catch
    errordlg('You must specify both model order and training points')
    stop %force it to stop here after throwing error report
end

if(isnan(order) || isnan(app_end))
    errordlg('You must specify both model order and training points')
    stop
end


ver_start=1;
try
t=handles.time;
catch
    errordlg('You must double click to download data');
end
ver_end = size(t,2);
if(app_end > ver_end)
  errordlg('Error: the training length exceeds the length of the data sequence');
end
Ts = t(2)-t(1);

beg = 1;
max_port=size(handles.y_temp_vec,1);
pt = floor(max_port/2);
offset = 293.15;%handles.y_temp_vec(1,1);

y_app = handles.y_temp_vec(:,app_start:app_end)-offset;
u_app = handles.u_p_vec(:,app_start:app_end);

y_ver = handles.y_temp_vec(:,ver_start:ver_end)-offset;
u_ver = handles.u_p_vec(:,ver_start:ver_end);

data_orig = iddata(y_app',u_app',Ts);
data_ver = iddata(y_ver',u_ver',Ts);
display('start subspace identification, please wait......');
h = waitbar(0,'Please wait...');

ChoiceOfSID = handles.current_choice
switch ChoiceOfSID
case 'N4SID'
    [ A B C D X0 SSM ] = subspace(u_app', y_app', t', order);
case 'PEM'
    [ A B C D X0 SSM ] = subspace_pem(u_app', y_app', t', order);
case 'SID_DET'
    [ A B C D X0 SSM ] = subspace_det(u_app', y_app', t', order);
case 'SUBID'
    [ A B C D X0 SSM ] = subspace_subid(u_app', y_app', t', order);
end
    

yTempSid(:,1:ver_end)=handles.y_temp_vec(:,1:ver_end)-offset;

%[yx,fit,xxo] = compare(data_ver,SSM);
tic
xh = X0;
for i=1:ver_end
%     if i==app_end
%         tic
%     end
  xh = A*xh + B*handles.u_p_vec(:,i);
  yh(:,i) = C*xh + D*handles.u_p_vec(:,i);
end
%yh=yx{1,1}.OutputData';
save('ABCDx0.mat', 'A','B','C','D','X0');
%handles.y_sid = yh{1,1}.OutputData(:,:);
handles.y_sid = yh;
simulation_time = toc

%save the temperature variable to the workspace
assignin('base','yTempSid',yTempSid);
assignin('base','t',t);
guidata(hObject, handles);

err_stats = zeros(max_port,2);
%length=size(yh{1,1}.OutputData,1);
length=size(yh,2);
for i = 1:max_port
    %err_diff = abs(yh{1,1}.OutputData(:,i) - y_ver(i,:)');
    err_diff = abs(yh(i,:) - y_ver(i,:));
    err_per = err_diff./(y_ver(i,:)+0*273.15+1);
    err_stats(i,1) = 100*mean(err_per);
    err_stats(i,2) = 100*std(err_per);
end
handles.err = max(err_stats(:,1));
%handles.err = 'transfer';
handles.identity=1;
handles.model_ss = SSM;
  Ac = SSMC.A;
  Bc = SSMC.B;
  Cc = SSMC.C;
  Dc = SSMC.D;
  Xc0 = SSMC.X0;
  
  [V,Lambda] = eig(Ac);
  AA = Lambda;
  BB = inv(V)*Bc;
  CC = Cc*V;
  DD = Dc;
  KK = zeros(size(Ac,1),size(Cc,1) );
  XX0 = inv(V)*Xc0;
  Ts = 0;
  
  %%%%%find its frequency domain response
  Nord=size(Ac,1);
  Ns = 1000;
  s=2*pi*j*logspace(-3,2,Ns); 
 Dk=zeros(Ns,Nord); 
  for m=1:Nord
    Dk(:,m)=1./(s-diag(AA(m)));
  end 
  %%%when input 1 is excited
  inNum = size(Bc,2);
  unitIn=zeros(inNum,1);
  unitIn(1)=1;
  Yfit=[];
 for k=1:Ns
     Yfit(:,k)=(CC*inv(s(k)*eye(size(AA)) - AA)*BB+DD)*unitIn;
 end
 handles.Yfit = Yfit;
 handles.s = s;


guidata(hObject, handles);
close(h);
display('Identification is finished.');




% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = get_var_names(handles);
 %figure; %figure(gcf)
 %axes(handles.SID_Result);

 evalin('base',['plot(',x,',',y,')'],'errordlg(lasterr,''Error generating plots'',''modal'')')

% title('Comparison of temperature output');
xlabel('time(sec)');
ylabel('temperature(Celcius)');
y=['temperature response ',y];
title(y)




function order_input_Callback(hObject, eventdata, handles)
order = str2double(get(handles.order_input,'String'));
handles.order = order;
guidata(hObject, handles);





function time_point_input_Callback(hObject, eventdata, handles)
% hObject    handle to time_point_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_point_input as text
%        str2double(get(hObject,'String')) returns contents of time_point_input as a double
app_end = str2double(get(handles.time_point_input,'String'));
handles.app_end = app_end;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function time_point_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_point_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in update_button.
function update_button_Callback(hObject, eventdata, handles)
% hObject    handle to update_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flag = handles.identity;
if(flag == 1 | flag ==2 )
    yTempSid = handles.y_sid+handles.y_temp_vec(1,1) - 273.15;
    fid = fopen('./matrixrow2vector.m','w');
    for i=1:size(yTempSid,1)
      fprintf(fid, 'y%d=yTempSid(%d,:)'';\n',i,i);
      fprintf(fid, 'assignin(\''base\'',\''y%d\'',y%d);\n',i,i);
%fprintf(fid, 'assignin(\'base','y_%d',y_%d);\n',i,i);
%fprintf(fid, 'u_%d=u_p_vec(%d,:);\n',i,i);
    end
    fclose(fid);
    run matrixrow2vector 
    %size(y1)
    update_listbox(handles)
    delete matrixrow2vector.*
else
update_listbox(handles)
end

function update_listbox(handles)
vars = evalin('base','who');
set(handles.listbox2,'String',vars)

function [var1,var2] = get_var_names(handles)
% Returns the names of the two variables to plot
list_entries = get(handles.listbox2,'String');
index_selected = get(handles.listbox2,'Value');
if length(index_selected) ~= 2
    errordlg('You must select two variables','Incorrect Selection','modal')
    
else
   try    
    var1 = list_entries{index_selected(1)};
    var2 = list_entries{index_selected(2)};
   catch
    errordlg('Variables are not unknown, please select the correct variables');
   end
end 
 





function Average_err_Callback(hObject, eventdata, handles)
% hObject    handle to Average_err (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Average_err as text
%        str2double(get(hObject,'String')) returns contents of Average_err as a double



% --- Executes during object creation, after setting all properties.
function Average_err_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Average_err (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in get_error.
function get_error_Callback(hObject, eventdata, handles)
flag = handles.identity;
if(flag==1 | flag==2)
set(handles.Average_err,'String',handles.err);
else
 errordlg('You must run SID first')
end
% hObject    handle to get_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in update_file.
function update_file_Callback(hObject, eventdata, handles)
% hObject    handle to update_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_file(pwd, handles)

function update_file(dir_path, handles)
cd (dir_path)
dir_struct = dir(dir_path);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.file_names = sorted_names;
handles.is_dir = [dir_struct.isdir];
handles.sorted_index = sorted_index;
guidata(handles.figure1,handles)
set(handles.listbox1,'String',handles.file_names,...
	'Value',1)





% --- Executes on button press in bodeplot_button.
function bodeplot_button_Callback(hObject, eventdata, handles)
flag = handles.identity;
if(flag==1)
sys = d2c(tf(handles.model_ss,'m')); 
try
m=handles.InPortNum;
n=handles.OutPortNum;
bode(sys(n,m),{1e-4,300});
catch
    errordlg('The input port or output port is not correctly specified');
end
elseif (flag==2)
    errordlg('Bode plot is not a proper representation for picewise linear model');
else
   errordlg('You must run SID first')
end



function In_port_Callback(hObject, eventdata, handles)
InPortNum = str2double(get(handles.In_port,'String'));
handles.InPortNum = InPortNum;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function In_port_CreateFcn(hObject, eventdata, handles)
% hObject    handle to In_port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Out_port_Callback(hObject, eventdata, handles)
OutPortNum = str2double(get(handles.Out_port,'String'));
handles.OutPortNum = OutPortNum;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Out_port_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Out_port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in print_tf_button.
function print_tf_button_Callback(hObject, eventdata, handles)
flag = handles.identity;
if(flag==1)
sys = d2c(tf(handles.model_ss,'m')); 
try
m=handles.InPortNum;
n=handles.OutPortNum;
sys(n,m)
catch
    errordlg('The input port or output port is not correctly specified');
end
elseif(flag==2)
     errordlg('Transfer function is not a proper representation for picewise linear model');
else
       
   errordlg('You must run SID first');
end





function selected_out_Callback(hObject, eventdata, handles)
% hObject    handle to selected_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of selected_out as text
%        str2double(get(hObject,'String')) returns contents of selected_out as a double
SelectOut = str2double(get(handles.selected_out,'String'));
handles.SelectOut = SelectOut;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function selected_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selected_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_comparison_button.
function plot_comparison_button_Callback(hObject, eventdata, handles)
flag = handles.identity;
if(flag==1 | flag==2)
% sys = d2c(tf(handles.model_ss,'m')); 
try
p=handles.SelectOut;
catch
    errordlg('Please specify the output port');
    stop;
end
 if(p>size(handles.y_temp_vec,1))
    errordlg('The output port is not correctly specified');
    stop;
 end
  try
   handles.y_temp_vec(p,:);
  catch
    errordlg('Please specify the output port');
    stop; 
  end
   figure; plot(handles.time, handles.y_temp_vec(p,:)-273.15,'g','linewidth',2); 
   hold;
   plot(handles.time, handles.y_sid(p,:)+handles.y_temp_vec(1,1)-273.15,'--r','linewidth',2);
   xlabel('time(sec)');
   ylabel('temperature(Celcius)');
   fid=fopen('title_fig.m','w');
   fprintf(fid, 'title(''temperature response at Output%d'')\n',p);
   fprintf(fid,'legend(''reference'',''sid'');\n');
   fclose(fid);
   run title_fig
   delete title_fig.m
   
else
   errordlg('You must run SID first')
end





% --- Executes on selection change in choose_method_popup.
function choose_method_popup_Callback(hObject, eventdata, handles)
% hObject    handle to choose_method_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns choose_method_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from choose_method_popup
val = get(hObject,'Value');
str = get(hObject, 'String');
%by default N4SID is in use
%handles.current_choice = 'N4SID';
switch str{val};
    case 'N4SID' % User selects peaks
	   handles.current_choice = 'N4SID';
    case 'PEM'
        handles.current_choice = 'PEM';
    case 'SID_DET'
        handles.current_choice = 'SID_DET';
    case 'SUBID'
        handles.current_choice = 'SUBID';
end
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function choose_method_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to choose_method_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function section_length_input_Callback(hObject, eventdata, handles)
% hObject    handle to section_length_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of section_length_input as text
%        str2double(get(hObject,'String')) returns contents of section_length_input as a double
SectLen = str2double(get(handles.section_length_input,'String'));
handles.SectLen = SectLen;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function section_length_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to section_length_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function shared_length_input_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of shared_length_input as text
%        str2double(get(hObject,'String')) returns contents of shared_length_input as a double
SharedLen = str2double(get(handles.shared_length_input,'String'));
handles.SharedLen = SharedLen;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function shared_length_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shared_length_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PWSID_button.
function PWSID_button_Callback(hObject, eventdata, handles)
% hObject    handle to PWSID_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
%load prbs_siso6_2.mat;
try
T_start=handles.Tstart;
T_step=handles.Tstep;
catch
    errordlg('You must specify T_start and T_step to partition the temperature into different ranges');
end
app_start = 1;
try
order = handles.order;
app_tot = handles.app_end;
catch
    errordlg('You must specify both model order and training points');
    stop
end

try
app_end = handles.SectLen;
app_offset = handles.SharedLen;
catch
    errordlg('You must specify section length and shared length for piecewise SID');
    stop
end
app_len = app_end%4800;
try
offset = 293.15; %handles.y_temp_vec(1,1);
catch
    errordlg('You must double click to download the data for piecewise SID');
    stop
end

ver_start=1;
ver_end = size(handles.y_temp_vec,2);
if(app_tot > ver_end)
   errordlg('Error: the training length exceeds the length of the data sequence')
   
end

if(mod(app_tot,app_len) | app_tot==app_len)
    errordlg('the total length should be integral(>1) of the section length');
    stop;
else
NumPw = app_tot/(app_len-app_offset);
end

t=handles.time;
Ts = t(2)-t(1);


display('start subspace identification, please wait......');
h = waitbar(0,'Please wait...');

y_app1 = handles.y_temp_vec(:,app_start:app_end)- offset;
u_app1 = handles.u_p_vec(:,app_start:app_end);
[ A1 B1 C1 D1 X01 SSM1 ] = subspace(u_app1', y_app1', t(app_start:app_end)', order);
xh1(:,1)=X01;
for i=1:app_len
  xh1(:,i+1) = A1*xh1(:,i) + B1*handles.u_p_vec(:,i);
  yh(:,i) = C1*xh1(:,i) + D1*handles.u_p_vec(:,i);
end
app_start = app_start + app_len - app_offset ;
app_end = app_start + app_len-1;
iter1=0;

T_pre = eye(size(A1));
xh_pre= xh1(:,end-app_offset:end-1);
for k=2:NumPw-1
y_app = handles.y_temp_vec(:,app_start:app_end)- offset;
u_app = handles.u_p_vec(:,app_start:app_end);
[ A(:,:,k) B(:,:,k) C(:,:,k) D(:,:,k) X0(:,k) SSM ] = subspace(u_app', y_app', t(app_start:app_end)', order);
xh(:,1)=X0(:,k);
for i=1:app_len
  xh(:,i+1) = reshape(A(:,:,k),size(A1))*xh(:,i) + reshape(B(:,:,k),size(B1))*handles.u_p_vec(:,i+app_start-1);
  yh(:,i+app_start-1) = reshape(C(:,:,k),size(C1))*xh(:,i) + reshape(D(:,:,k),size(D1))*handles.u_p_vec(:,i+app_start-1);
end

T = xh(:,1:app_offset)/xh_pre*T_pre;
T_inverse = inv(T);

At(:,:,k) = T_inverse*reshape(A(:,:,k), size(A1))*T;
Bt(:,:,k) = T_inverse*reshape(B(:,:,k),size(B1));
Ct(:,:,k) = reshape(C(:,:,k),size(C1))*T;
Dt(:,:,k) = reshape(D(:,:,k),size(D1));
T_pre = T;
xh_pre=xh(:,end-app_offset:end-1);
if(k<NumPw-1)
app_start = app_start + app_len - app_offset;
app_end = app_start + app_len-1;
end

end

%check the remaining data  and decide whether or not to build another model
if(app_end < app_tot)
    k=k+1;
    app_start = app_start + app_len - app_offset;
    y_app = handles.y_temp_vec(:,app_start:app_tot) - offset;
    u_app = handles.u_p_vec(:,app_start:app_tot);
    [ A(:,:,k) B(:,:,k) C(:,:,k) D(:,:,k) X0(:,k) SSM ] = subspace(u_app', y_app', t(app_start:app_tot)', order);
    xh(:,1)=X0(:,k);

    for i=1:app_tot-app_end+app_offset
    xh(:,i+1) = reshape(A(:,:,k),size(A1))*xh(:,i) + reshape(B(:,:,k),size(B1))*handles.u_p_vec(:,i+app_start-1);
    yh(:,i+app_start-1) = reshape(C(:,:,k),size(C1))*xh(:,i) + reshape(D(:,:,k),size(D1))*handles.u_p_vec(:,i+app_start-1);
    end
    T = xh(:,1:app_offset)/xh_pre*T_pre;
    T_inverse = inv(T);

    At(:,:,k) = T_inverse*reshape(A(:,:,k), size(A1))*T;
    Bt(:,:,k) = T_inverse*reshape(B(:,:,k),size(B1));
    Ct(:,:,k) = reshape(C(:,:,k),size(C1))*T;
    Dt(:,:,k) = reshape(D(:,:,k),size(D1));
  
end
%end
 save('piecewiseModel.mat','At','Bt','Ct','Dt','X01'); 
app_end=app_tot;


y_ver = handles.y_temp_vec(:,app_start:ver_end) - offset;
u_ver = handles.u_p_vec(:,app_start:ver_end)';


 xh=zeros(order,ver_end);
 xh(:,1) = X01;
 portN = size(handles.y_temp_vec,1);
 yh=zeros(portN,ver_end);
 yh(:,1) = C1*xh(:,1) + D1*handles.u_p_vec(:,1);
 %iter1=0; iter2=0; iter3=0; iter4=0; iter5=0; iter6=0; iter7=0; iter8=0; iter9=0; iter10=0; iter11=0; 
 %setup the simulation data partition

 fid2=fopen('sim_model.m','w');
  fprintf(fid2,'function [yh]=sim_model(T_start,T_step,t_0,t_end,app_end, offset,  yh, xh, At, Bt, Ct, Dt, p_vec, A1, B1, C1, D1)\n');
 for m=2:k
   fprintf(fid2, 'A%dt=reshape(At(:,:,%d),size(A1)); B%dt=reshape(Bt(:,:,%d),size(B1)); C%dt=reshape(Ct(:,:,%d),size(C1)); D%dt=reshape(Dt(:,:,%d),size(D1)); \n',m,m,m,m,m,m,m,m);
 end
 fprintf(fid2,'for i=t_0:t_end\n');
 %fprintf(fid2,'yh_m = y_temp_vec(1,i)-offset+20;\n');
 fprintf(fid2,'if(i==1) yh_m = sum(yh(:,1))/25+20;\n');
 fprintf(fid2,'else yh_m=sum(yh(:,i-1))/25+20; end\n');
 fprintf(fid2, 'if(yh_m<T_start)\n');
 fprintf(fid2,'xh(:,i+1) = A1*xh(:,i) + B1*p_vec(:,i);\n');
 fprintf(fid2,'yh(:,i) = C1*xh(:,i) + D1*p_vec(:,i);\n');
 
 for m=2:k %k is the number of the model
     if(m<k)
      fprintf(fid2, 'elseif(yh_m>T_start+T_step*(%d-2) && yh_m< T_start+T_step*(%d-1) )\n',m,m);
     else
      fprintf(fid2, 'elseif(yh_m> T_start+T_step*(%d-2) )\n',m);  
     end     
      fprintf(fid2,'xh(:,i+1) = A%dt*xh(:,i) + B%dt*p_vec(:,i);\n',m,m);
      fprintf(fid2,'yh(:,i) = C%dt*xh(:,i) + D%dt*p_vec(:,i);\n',m,m);
 end
 fprintf(fid2,'end\n');
 fprintf(fid2,'end\n');
 fprintf(fid2,'yhs=yh;');
 fclose(fid2);
 tic
 t0=1;
 t_end= ver_end;

 yhs=sim_model(T_start, T_step, t0, t_end, app_end, offset,  yh, xh, At, Bt, Ct, Dt, handles.u_p_vec, A1, B1, C1, D1);
 simulation_time=toc

  err_diff = abs(handles.y_temp_vec(:,app_end:ver_end)-yhs(:,app_end:ver_end) - offset)./(handles.y_temp_vec(:,app_end:ver_end)-273.15);%-273.15);
  mean_err = 100*mean(err_diff');
  max_err = max(mean_err);

%plot(t(1:ver_end)',yhs(1,1:ver_end)+20,'--r','linewidth',2); hold on;

handles.y_sid = yhs;
handles.err=max_err;
handles.identity=2;
guidata(hObject, handles);
close(h);
display('Identification is finished.');


function Temp_init_input_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of Temp_init_input as text
%        str2double(get(hObject,'String')) returns contents of Temp_init_input as a double

Tstart = str2double(get(handles.Temp_init_input,'String'));
handles.Tstart = Tstart;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Temp_init_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Temp_init_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Temp_step_input_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of Temp_step_input as text
%        str2double(get(hObject,'String')) returns contents of Temp_step_input as a double
Tstep = str2double(get(handles.Temp_step_input,'String'));
handles.Tstep = Tstep;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Temp_step_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Temp_step_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function netlistName_Callback(hObject, eventdata, handles)
% hObject    handle to netlistName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of netlistName as text
%        str2double(get(hObject,'String')) returns contents of netlistName as a double
netName = get(handles.netlistName,'String');
% fid1=fopen(netName,'w');
% fclose(fid1);
handles.genNet = 0;
handles.netName = netName;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function netlistName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to netlistName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GenNetlist_button.
function GenNetlist_button_Callback(hObject, eventdata, handles)
% hObject    handle to GenNetlist_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  flag = handles.identity;
  if(flag == 1)
     SSMC = d2c(handles.model_ss);
  elseif (flag == 2)
      errordlg('cannot implement piecewise linear model');
      stop
  else
      errordlg('Run SID first!');
  end
  Ac = SSMC.A;
  Bc = SSMC.B;
  Cc = SSMC.C;
  Dc = SSMC.D;
  Xc0 = SSMC.X0;
  
  [V,Lambda] = eig(Ac);
  AA = Lambda;
  BB = inv(V)*Bc;
  CC = Cc*V;
  DD = Dc;
  KK = zeros(size(Ac,1),size(Cc,1) );
  XX0 = inv(V)*Xc0;
  Ts = 0;
  
  %%%%%find its frequency domain response
  Nord=size(Ac,1);
  Ns = 1000;
  s=2*pi*j*logspace(-3,2,Ns); 
 Dk=zeros(Ns,Nord); 
  for m=1:Nord
    Dk(:,m)=1./(s-diag(AA(m)));
  end 
  %%%when input 1 is excited
  inNum = size(Bc,2);
  unitIn=zeros(inNum,1);
  unitIn(1)=1;
  Yfit=[];
 for k=1:Ns
     Yfit(:,k)=(CC*inv(s(k)*eye(size(AA)) - AA)*BB+DD)*unitIn;
 end
 handles.Yfit = Yfit;
 handles.s = s;
 handles.genNet = 1;
 guidata(hObject,handles);
 
cindex=zeros(1,Nord);
for m=1:Nord 
  if imag(Lambda(m,m))~=0  
    if m==1 
      cindex(m)=1;
    else
      if cindex(m-1)==0 || cindex(m-1)==2
        cindex(m)=1; cindex(m+1)=2; 
      else
        cindex(m)=2;
      end
    end 
  end
end

  n=0;
  for m=1:Nord
    n=n+1;
    if cindex(m)==1
      a=AA(n,n); a1=real(a); a2=imag(a);
      c=CC(:,n); c1=real(c); c2=imag(c);
      b=BB(n,:); b1=2*real(b); b2=-2*imag(b);
      Ablock=[a1 a2;-a2 a1];
   
      AA(n:n+1,n:n+1)=Ablock;
      CC(:,n)=c1;
      CC(:,n+1)=c2;
      BB(n,:)=b1;
      BB(n+1,:)=b2;
    end
  end
disp('The system matrix A is converted to complex');
 
 sim_opts.scanType = 'DEC';
 sim_opts.scanNum = '100';
 sim_opts.freqStart = num2str(1e-3);
 sim_opts.freqEnd = num2str(100); 
 circuit_name = strtok(handles.netName,'.');
 [Gi, Ci, Bi, Ui] = rlz_vf3(AA, BB, CC, DD, size(BB,2), size(CC,1), handles.netName, circuit_name, 1, 0, 1, sim_opts);
  
  disp('File write to spice netlist');
   
   

% --- Executes on button press in RunSpice_button.
function RunSpice_button_Callback(hObject, eventdata, handles)
% hObject    handle to RunSpice_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
     disp('Start Spice simulatin...');
    % [sx,wx] = dos( ['hspice ckt.cir -o cir.lis'] )
    fid_cir = fopen('RunHspice.m','w');
    fprintf(fid_cir, '[sx,wx] = dos([''hspice %s -o cir.lis''])',handles.netName);
    fclose(fid_cir);
    run RunHspice;


% --- Executes on button press in plotAC_button.
function plotAC_button_Callback(hObject, eventdata, handles)
% hObject    handle to plotAC_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
     x=loadsig('cir.ac0');
     fr=evalsig(x,'HERTZ');
     PortN = handles.PortNum;
     PortS = num2str(handles.PortNum);
     try
     q=evalsig(x,strcat('o_port_',PortS));
     catch
         errordlg('The port number exceed the range');
         stop
     end
     figure; 
     if(handles.genNet==1)
     h4 = semilogx((fr),db(abs(q)),'g');
     hold on;
     s = handles.s;
     Yfit = handles.Yfit;
     h5 = semilogx(abs(s/(2*pi*j)),db(abs(Yfit(PortN,:))),'--r');
     legend([h4(1) h5(1)],'Spice Simulation','N4SID');
     xlabel('Frequency [Hz]')
     ylabel('Gain [db]');
     title('Transfer functions');
     else
     h4 = semilogx((fr),db(abs(q)),'g');
     legend([h4(1) ],'Spice Simulation');
     xlabel('Frequency [Hz]')
     ylabel('Gain [db]');
     title('Transfer functions'); 
     end


function OutPort_Callback(hObject, eventdata, handles)
% hObject    handle to OutPort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PortNum = str2double(get(handles.OutPort,'String'));
handles.PortNum = PortNum;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function OutPort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutPort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotTran_button.
function plotTran_button_Callback(hObject, eventdata, handles)
% hObject    handle to plotTran_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
     x=loadsig('cir.tr0');
     t=evalsig(x,'TIME');
     PortS = num2str(handles.PortNum);
     try
     q=evalsig(x,strcat('o_port_',PortS));
     catch
         errordlg('The port number exceed the range');
         stop
     end
     plot(t,q)
     xlabel('time(sec)')
     ylabel('temp(Celcius)');
     legend('Spice Simulation');
     

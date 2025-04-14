function varargout = Global_Apparent_Density_Mapping(varargin)
% GLOBAL_APPARENT_DENSITY_MAPPING MATLAB code for Global_Apparent_Density_Mapping.fig
%      GLOBAL_APPARENT_DENSITY_MAPPING, by itself, creates a new GLOBAL_APPARENT_DENSITY_MAPPING or raises the existing
%      singleton*.
%
%      H = GLOBAL_APPARENT_DENSITY_MAPPING returns the handle to a new GLOBAL_APPARENT_DENSITY_MAPPING or the handle to
%      the existing singleton*.
%
%      GLOBAL_APPARENT_DENSITY_MAPPING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GLOBAL_APPARENT_DENSITY_MAPPING.M with the given input arguments.
%
%      GLOBAL_APPARENT_DENSITY_MAPPING('Property','Value',...) creates a new GLOBAL_APPARENT_DENSITY_MAPPING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Global_Apparent_Density_Mapping_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Global_Apparent_Density_Mapping_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Global_Apparent_Density_Mapping

% Last Modified by GUIDE v2.5 02-Apr-2025 19:34:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Global_Apparent_Density_Mapping_OpeningFcn, ...
                   'gui_OutputFcn',  @Global_Apparent_Density_Mapping_OutputFcn, ...
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


% --- Executes just before Global_Apparent_Density_Mapping is made visible.
function Global_Apparent_Density_Mapping_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Global_Apparent_Density_Mapping (see VARARGIN)

% Choose default command line output for Global_Apparent_Density_Mapping
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Global_Apparent_Density_Mapping wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Global_Apparent_Density_Mapping_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global upboundry lowboundry gravity  accumulated_density interval lon_min lon_max lat_min lat_max

lon_min=str2double(get(handles.edit5,'string'));
lon_interval=str2double(get(handles.edit6,'string'));
lon_max=str2double(get(handles.edit7,'string'));
lat_min=str2double(get(handles.edit8,'string'));
lat_interval=str2double(get(handles.edit9,'string'));
lat_max=str2double(get(handles.edit10,'string'));
nlon=(lon_max-lon_min)/lon_interval+1;
nlat=(lat_max-lat_min)/lat_interval+1;
Lat=linspace(lat_max,lat_min,nlat); Lon=linspace(lon_min,lon_max,nlon);
interval=lon_interval;
initial_density=0;
height=str2double(get(handles.edit2,'string'));
max_iter= str2double(get(handles.edit3,'string'));  
tol=str2double(get(handles.edit4,'string'));  
alpha =str2double(get(handles.edit1,'string'));  
nmax=str2double(get(handles.edit12,'string'));

up_boundary=upboundry;low_boundary=lowboundry;

 

if handles.radiobutton1.Value == 1
  gravity_attraction=gravity;  
[accumulated_density,gravity_rms_values] = accumulation_density_by_attraction(up_boundary, low_boundary, height, gravity_attraction, initial_density, max_iter, tol, alpha,Lon,Lat,nmax);
else 
  gravity_gradient=gravity;  
[accumulated_density,gravity_rms_values] = accumulation_density_by_gradient(up_boundary, low_boundary, height, gravity_gradient, initial_density, max_iter, tol, alpha,Lon,Lat,nmax);
end
fmin=min(accumulated_density(:,max_iter+2));
fmax=max(accumulated_density(:,max_iter+2));
fmean=mean(accumulated_density(:,max_iter+2));
fstd=std(accumulated_density(:,max_iter+2));
stastic=[fmin,fmax,fmean,fstd];
set(handles.uitable5,'data',stastic);

fresult(:,1:2)=accumulated_density(:,1:2);
fresult(:,3)=accumulated_density(:,max_iter+2);
set(handles.uitable4,'data',fresult);    

axes(handles.axes1);
plot(1:max_iter, gravity_rms_values, '-o');
xlabel('Iteration', 'FontName', 'Times New Roman');
ylabel('RMS Error', 'FontName', 'Times New Roman');
if handles.radiobutton1.Value == 1
title('Gravity Anomaly RMS Curve', 'FontName', 'Times New Roman');
else
title('Gravity Gradient RMS Curve', 'FontName', 'Times New Roman');
end
grid on;


axes(handles.axes2);
scatter(accumulated_density(:,1),accumulated_density(:,2),1,accumulated_density(:,max_iter+2));
cb=colorbar;colormap(jet); 
set(get(cb, 'Title'), 'String', 'g/cm^3', 'FontName', 'Times New Roman');
title(['Density Image of Iteration ', num2str(max_iter)], 'FontName', 'Times New Roman');

xlim([lon_min lon_max]);
ylim([lat_min lat_max]);
xlabel('Longitude (degree)', 'FontName', 'Times New Roman');
ylabel('Latitude (degree)', 'FontName', 'Times New Roman');







% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


resultdata=get(handles.uitable4,'data');
[filename,pathname,c]=uiputfile('*.txt','save');
if c==1
file=[pathname,filename];
dlmwrite(file,resultdata);
helpdlg('Successfully saved!')
end

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitable4,'data','');
set(handles.uitable5,'data','');


function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global accumulated_density  lon_min lon_max lat_min lat_max
N=str2double(get(handles.edit11,'string'));
set(handles.uitable4,'data','');
set(handles.uitable5,'data','');
cla(handles.axes2, 'reset');

fmin=min(accumulated_density(:,N+2));
fmax=max(accumulated_density(:,N+2));
fmean=mean(accumulated_density(:,N+2));
fstd=std(accumulated_density(:,N+2));
stastic=[fmin,fmax,fmean,fstd];
set(handles.uitable5,'data',stastic);

fresult(:,1:2)=accumulated_density(:,1:2);
fresult(:,3)=accumulated_density(:,N+2);
set(handles.uitable4,'data',fresult);    

axes(handles.axes2);
scatter(accumulated_density(:,1),accumulated_density(:,2),1,accumulated_density(:,N+2));
cb=colorbar;colormap(jet); 
set(get(cb, 'Title'), 'String', 'g/cm^3', 'FontName', 'Times New Roman');
title(['Density Image of Iteration ', num2str(N)], 'FontName', 'Times New Roman');
xlim([lon_min lon_max]);
ylim([lat_min lat_max]);
% xticks(lon_min:60:lon_max);
% yticks(lat_min:30:lat_max);
xlabel('Longitude (degree)', 'FontName', 'Times New Roman');
ylabel('Latitude (degree)', 'FontName', 'Times New Roman');

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global interval
h = figure('visible','off');
new_axes = copyobj(handles.axes2,h);
set(new_axes,'Units','default','Position','default');
hScatter = findobj(new_axes, 'Type', 'Scatter');
    p=2*interval;
    set(hScatter, 'SizeData',p );
    cb=colorbar;
    set(get(cb,'title'),'string','g/cm^3');
    colormap(jet)
    [filename, pathname] = uiputfile({'*.png';'*.bmp';'*.jpg'},'The picture is saved as');
if filename ~= 0
        file = strcat(pathname,filename);
        saveas(h,file);
        msgbox('The image has been successfully saved','Attention','help');
else
        msgbox('The operation has been canceled','Attention','warn');
end

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes2, 'reset');

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_f_handle=figure('visible','off');
new_axes=copyobj(handles.axes1,new_f_handle);
set(new_axes,'units','default','Position','default');
[filename,pathname]=uiputfile({'*.png'},'save picture as');
if ~filename
    return
else
    file=strcat(pathname,filename);
    print(new_f_handle,'-djpeg',file);
end
delete(new_f_handle);
h=dialog('Name','Save data','Position',[200 200 200 70]);
uicontrol('Style','text','Units','pixels','Position',[50 40 120 20],'FontSize',10,...
    'Parent',h,'String','save done'); 
uicontrol('Units','pixels','Position',[80 10 50 20],'FontSize',10,...
    'Parent',h,'String','OK','Callback','delete(gcf)');


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1, 'reset');

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gravity
[filename,pathname]=uigetfile('*.txt','open file');
if filename==0
    return
else
    file=[pathname,filename];
    data=csvread(file);
    set(handles.uitable3,'data',data);
end 
gravity=data;
int=data(2,1)-data(1,1);
nmax=180/(int);
 set(handles.edit12, 'string', nmax)

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitable3,'data','');

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global upboundry lowboundry
[filename,pathname]=uigetfile('*.txt','open file');
if filename==0
    return
else
    file=[pathname,filename];
    data=csvread(file);

table_data = get(handles.uitable1,'data');

if size(table_data,1)==4||size(table_data,1)==0
    table_data=[];
end
table_data(:,[1,2])=data(:,[1,2]);
table_data=[table_data,data(:,3)];
set(handles.uitable1,'data',table_data);
end 
if size(table_data,2)==4
upboundry(:,1:2)=table_data(:,1:2);
upboundry(:,3)=table_data(:,3);
lowboundry(:,1:2)=table_data(:,1:2);
lowboundry(:,3)=table_data(:,4);
end
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitable1,'data','');



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

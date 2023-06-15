function varargout = Multilayer_PlaneWave(varargin)
% MULTILAYER_PLANEWAVE MATLAB code for Multilayer_PlaneWave.fig
%      MULTILAYER_PLANEWAVE, by itself, creates a new MULTILAYER_PLANEWAVE or raises the existing
%      singleton*.
%
%      H = MULTILAYER_PLANEWAVE returns the handle to a new MULTILAYER_PLANEWAVE or the handle to
%      the existing singleton*.
%
%      MULTILAYER_PLANEWAVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTILAYER_PLANEWAVE.M with the given input arguments.
%
%      MULTILAYER_PLANEWAVE('Property','Value',...) creates a new MULTILAYER_PLANEWAVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Multilayer_PlaneWave_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Multilayer_PlaneWave_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Multilayer_PlaneWave

% Last Modified by GUIDE v2.5 04-Mar-2018 09:59:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Multilayer_PlaneWave_OpeningFcn, ...
                   'gui_OutputFcn',  @Multilayer_PlaneWave_OutputFcn, ...
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


% --- Executes just before Multilayer_PlaneWave is made visible.
function Multilayer_PlaneWave_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Multilayer_PlaneWave (see VARARGIN)

% Choose default command line output for Multilayer_PlaneWave
handles.output = hObject;

set(handles.edit_theta_sweep, 'enable', 'off');
set(handles.edit_layer_number2, 'enable', 'off');
set(handles.popupmenu_rt2, 'enable', 'off');
        
set(handles.edit_freq_sweep, 'enable', 'on');
set(handles.edit_layer_number, 'enable', 'on');
set(handles.popupmenu_rt, 'enable', 'on');
%%
eps_init_val = eval(get(handles.edit_eps_r, 'string'));
mu_init_val = eval(get(handles.edit_eps_r, 'string'));
freq_init_val = eval(get(handles.edit_freq, 'string'));

thickness_init_val = eval(get(handles.edit_thickness, 'string'));
x = eval(get(handles.edit_x, 'string'));
y = eval(get(handles.edit_y, 'string'));
%%
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Multilayer_PlaneWave wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Multilayer_PlaneWave_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_CST.
function pushbutton_CST_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_CST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

thickness = eval(get(handles.edit_thickness, 'string'))*1e-3;
eps_r = eval(get(handles.edit_eps_r, 'string'));
mu_r = eval(get(handles.edit_mu_r, 'string'));

freq = eval(get(handles.edit_freq_sweep, 'string'));
angle0 = str2double(get(handles.edit_theta, 'string'));

nLayer = length(thickness);

fmin = freq(1);
fmax = freq(end);

cst01 = actxserver('CSTStudio.application');
mws = cst01.invoke('newMWS');


%Set Solver Type

mws.invoke('ChangeSolverType', 'HF Frequency Domain');

FDsolver = mws.invoke('FDSolver');
FDsolver.invoke('AddToExcitationList', 'Zmax', 'TE(0,0);TM(0,0)');

%Set Mesh Type
mesh01 = mws.invoke('Mesh');
mesh01.invoke('MeshType', 'PBA');
mesh01.invoke('SetCreator', 'High Frequency');

% meshSetting01 = mesh01.invoke('MeshSetting');
% meshSetting01.invoke('SetMeshType', 'Tet');



%Set Uinits

units01 = mws.invoke('Units'); %Create object Units

units01.invoke('Geometry', 'mm');
units01.invoke('Frequency', 'GHz');
units01.invoke('Voltage', 'V');
units01.invoke('Time', 'ns');

%Set Background Material
bckgnd01 = mws.invoke('Background');

bckgnd01.invoke('Type', 'Normal');

%Set Boundary Condition
boundary01 = mws.invoke('Boundary');

boundary01.invoke('Xmin', 'unit cell');
boundary01.invoke('Xmax', 'unit cell');
boundary01.invoke('Ymin', 'unit cell');
boundary01.invoke('Ymax', 'unit cell');
boundary01.invoke('Zmin', 'open');
boundary01.invoke('Zmax', 'open');
boundary01.invoke('SetPeriodicBoundaryAngles', angle0, '0');




%Set Frequency Range
solver01 = mws.invoke('Solver');

solver01.invoke('FrequencyRange', num2str(fmin), num2str(fmax));

%Create Model
build_structure(nLayer, eps_r, mu_r, thickness, 0, mws );


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sweep_type = get(handles.popupmenu_sweep, 'value');

switch sweep_type
    case 1 %Frequency Sweep
        freq = eval(get(handles.edit_freq_sweep, 'string'))*1e9;
        
        fmin = freq(1);
        fmax = freq(end);
        nf = length(freq);
        
        thickness = eval(get(handles.edit_thickness, 'string'))*1e-3;
        eps_r = eval(get(handles.edit_eps_r, 'string'));
        mu_r = eval(get(handles.edit_mu_r, 'string'));
        
        ang = str2double(get(handles.edit_theta, 'string'));
        
        nnLayer = length(thickness) ;
        
        eps_0 = 8.854187e-12;
        mu_0 = 4*pi*1e-7;
        
        k0 = 2*pi*freq*sqrt(eps_0 * mu_0);
        
         lastLayer_Val = get(handles.popupmenu_bc, 'value');

        switch lastLayer_Val
            case 1
                y_end = 0;
            case 2
                y_end = -1;
            case 3
                y_end = 1;
            case 4
                y_end = 0;
        end
        
        for jj = 1:length(freq)
            for ii = 1:nnLayer
                
                 %Calculate kd
                d(ii) = sum(thickness(1:ii));
    
                if (eps_r(ii)<0 && mu_r(ii)<0)
                    kd(ii, jj) = -k0(jj) * sqrt(abs(eps_r(ii)) * abs(mu_r(ii)));
                else
                    kd(ii, jj) = k0(jj) * sqrt(eps_r(ii) * mu_r(ii));
                end
                
               if ii == 1
                        theta(ii, jj) = ang;
               else
                        theta(ii, jj) = asind(kd(ii-1, jj) / kd(ii, jj) * sind(theta(ii - 1, jj)));
               end
                
            end
            
            
            for ii = nnLayer:-1:2
                if ii == nnLayer
        
                    y(ii - 1, jj) = recursive_model(y_end, d(ii - 1), kd(ii-1, jj)*cosd(theta(ii - 1, jj)), kd(ii, jj)*cosd(theta(ii, jj)), mu_0*mu_r(ii-1), mu_0*mu_r(ii));
        
                else 
                    y(ii - 1, jj) = recursive_model(y(ii, jj), d(ii - 1), kd(ii-1, jj)*cosd(theta(ii - 1, jj)), kd(ii, jj)*cosd(theta(ii, jj)), mu_0*mu_r(ii-1), mu_0*mu_r(ii));
                end
            end
            
            for ii = 1:nnLayer
    
                if ii == 1
                    mag_plus(ii, jj) = 1;
                    mag_minus(ii, jj) = mag_plus(ii, jj)*y(ii, jj);
                else
                    [mag_plus(ii, jj), mag_minus(ii, jj)] = forward_calculation(mag_plus(ii - 1, jj), mag_minus(ii - 1, jj), d(ii - 1), kd(ii-1, jj)*cosd(theta(ii - 1, jj)), kd(ii, jj)*cosd(theta(ii, jj)), mu_0*mu_r(ii-1), mu_0*mu_r(ii));
                end
    
            end  
            
            
            
        end
        
        selectRT = get(handles.popupmenu_rt, 'value');
        
        switch selectRT
                case 1
                    R = mag_minus(1, 1:end)./mag_plus(1,1:end);
                    T = mag_plus(nnLayer, 1:end) ./ mag_plus(1,1:end);
            
                    figure;
                    subplot(2,1,1);            
%                     plot(freq/1e9, 10*log10(abs(R).^2), 'LineWidth', 3, 'color', 'red');
                    plot(freq/1e9, abs(R), 'LineWidth', 3, 'color', 'red');
                    %axis([fmin/1e9 fmax/1e9 0 1]);
                    xlabel('Frequency (GHz)');
                    ylabel('Reflection Coefficient');
                    ylim([0 1]);
                    grid on;
            
                    subplot(2,1,2);
                    plot(freq/1e9, abs(T), 'LineWidth', 3, 'color', 'blue');
                    %plot(freq/1e9, 10*log10(abs(T).^2), 'LineWidth', 3, 'color', 'blue');
                    %axis([fmin/1e9 fmax/1e9 0 2]);
                    xlabel('Frequency (GHz)');
                    ylabel('Transmission Coefficient');
                    ylim([0 1]);
                    grid on;
            
                case 2
                    set(handles.edit_layer_number, 'enable', 'on');
                    layerNumber = str2double(get(handles.edit_layer_number, 'string'));
            
                % Plot Magnitude of E_plus
                    figure;
            
            
                    subplot(2,1,1);            
                    plot(freq/1e9, abs(mag_plus(layerNumber, 1:end)), 'LineWidth', 3, 'color', 'red'); %!!!!!!!!!!!!! Write Exception Handling Routine
                    
                    abs(mag_plus(layerNumber, 1:end))
                    axis([fmin/1e9 fmax/1e9 -2 2]);
            
                    strLabel = strcat('| E^{+}_{layer',num2str(layerNumber), '} |');
                    xlabel('Frequency (GHz)');
                    ylabel(strLabel);
                    grid on;
            
                    % Plot Magnitude of E_minus
                    subplot(2,1,2);
                    plot(freq/1e9, abs(mag_minus(layerNumber, 1:end)), 'LineWidth', 3, 'color', 'blue'); %!!!!!!!!!!!!!!!!!! Write Exception Handling Routine
                    axis([fmin/1e9 fmax/1e9 -2 2]);
            
                    strLabel = strcat('| E^{-}_{layer',num2str(layerNumber), '} |');
                    xlabel('Frequency (GHz)')
                    ylabel(strLabel);
                    grid on;
            
         end
        
        
        
    case 2 %Angular Sweep
        
        thickness = eval(get(handles.edit_thickness, 'string'))*1e-3;
        eps_r = eval(get(handles.edit_eps_r, 'string'));
        mu_r = eval(get(handles.edit_mu_r, 'string'));
        
        eps_0 = 8.854187e-12;
        mu_0 = 4*pi*1e-7;
        
        nnLayer = length(thickness) ;
        
        ang = eval(get(handles.edit_theta_sweep, 'string'));
        minTheta = ang(1);
        maxTheta = ang(end);
        nTheta = length(ang);
        
        kd = zeros(nnLayer, nTheta);
        
        freq = str2double(get(handles.edit_freq, 'string'))*1e9;
        
        k0 = 2*pi*freq*sqrt(eps_0 * mu_0);
        
        
        lastLayer_Val = get(handles.popupmenu_bc, 'value');

        switch lastLayer_Val
            case 1
                y_end = 0;
            case 2
                y_end = -1;
            case 3
                y_end = 1;
            case 4
                y_end = 0;
        end
        
        
         for jj = 1:length(ang)
    
            for ii = 1:nnLayer % Calculate parameters
        
                %Calculate kd
                d(ii) = sum(thickness(1:ii));
    
                if (eps_r(ii)<0 && mu_r(ii)<0)
                    kd(ii) = -k0 * sqrt(abs(eps_r(ii)) * abs(mu_r(ii)));
                else
                    kd(ii) = k0 * sqrt(eps_r(ii) * mu_r(ii));
                end
        
                %Phase Matching
    
                     if ii == 1
                        theta(ii, jj) = ang(jj);
                    else
                        theta(ii, jj) = asind(kd(ii-1) / kd(ii) * sind(theta(ii - 1, jj)));
                     end
             
            end
            
            %% Calculate Coefficient using Recursive Method
    
            for ii = nnLayer:-1:2
                if ii == nnLayer
        
                    y(ii - 1, jj) = recursive_model(y_end, d(ii - 1), kd(ii-1)*cosd(theta(ii - 1, jj)), kd(ii)*cosd(theta(ii, jj)), mu_0*mu_r(ii-1), mu_0*mu_r(ii));
        
                else 
                    y(ii - 1, jj) = recursive_model(y(ii, jj), d(ii - 1), kd(ii-1)*cosd(theta(ii - 1, jj)), kd(ii)*cosd(theta(ii, jj)), mu_0*mu_r(ii-1), mu_0*mu_r(ii));
                end
            end

            for ii = 1:nnLayer
    
                if ii == 1
                    mag_plus(ii, jj) = 1;
                    mag_minus(ii, jj) = mag_plus(ii, jj)*y(ii, jj);
                else
                    [mag_plus(ii, jj), mag_minus(ii, jj)] = forward_calculation(mag_plus(ii - 1, jj), mag_minus(ii - 1, jj), d(ii - 1), kd(ii-1)*cosd(theta(ii - 1, jj)), kd(ii)*cosd(theta(ii, jj)), mu_0*mu_r(ii-1), mu_0*mu_r(ii));
                end
    
            end  
        
        
        end
        
        selectRT = get(handles.popupmenu_rt2, 'value');
        
         switch selectRT
                case 1
                    R = mag_minus(1, 1:end)./mag_plus(1,1:end);
                    T = mag_plus(nnLayer, 1:end) ./ mag_plus(1,1:end);
            
                    figure;
                    subplot(2,1,1);            
                    plot(ang, abs(R), 'LineWidth', 3, 'color', 'red');
                    axis([minTheta maxTheta 0 5]);
                    xlabel('Incident Angle (deg)');
                    ylabel('Reflection Coefficient');
                    ylim([0 1]);
                    grid on;
            
                    subplot(2,1,2);
                    plot(ang, abs(T), 'LineWidth', 3, 'color', 'blue');
                    axis([minTheta maxTheta 0 5]);
                    xlabel('Incident Angle (deg)');
                    ylabel('Transmission Coefficient');
                    ylim([0 1]);
                    grid on;
            
                case 2
                    set(handles.edit_layer_number2, 'enable', 'on');
                    layerNumber = str2double(get(handles.edit_layer_number2, 'string'));
            
            % Plot Magnitude of E_plus
                    figure;
            
            
                    subplot(2,1,1);            
                    plot(ang, abs(mag_plus(layerNumber, 1:end)), 'LineWidth', 3, 'color', 'red'); %!!!!!!!!!!!!! Write Exception Handling Routine
                    axis([minTheta maxTheta -2 2]);
            
                    strLabel = strcat('| E^{+}_{layer',num2str(layerNumber), '} |');
                    xlabel('Incident Angle (deg)');
                    ylabel(strLabel);
                    grid on;
            
                    % Plot Magnitude of E_minus
                    subplot(2,1,2);
                    plot(ang, abs(mag_minus(layerNumber, 1:end)), 'LineWidth', 3, 'color', 'blue'); %!!!!!!!!!!!!!!!!!! Write Exception Handling Routine
                    axis([minTheta maxTheta -2 2]);
            
                    strLabel = strcat('| E^{-}_{layer',num2str(layerNumber), '} |');
                    xlabel('Incident Angle (deg)')
                    ylabel(strLabel);
                    grid on;
            
          end
        
        
        
        
end


% --- Executes on selection change in popupmenu_sweep.
function popupmenu_sweep_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_sweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sweep_type = get(handles.popupmenu_sweep, 'value');
switch sweep_type
    case 1
        set(handles.edit_theta_sweep, 'enable', 'off');
        set(handles.edit_layer_number2, 'enable', 'off');
        set(handles.popupmenu_rt2, 'enable', 'off');
        
        set(handles.edit_freq_sweep, 'enable', 'on');
        set(handles.edit_layer_number, 'enable', 'on');
        set(handles.popupmenu_rt, 'enable', 'on');
    case 2
        set(handles.edit_theta_sweep, 'enable', 'on');
        set(handles.edit_layer_number2, 'enable', 'on');
        set(handles.popupmenu_rt2, 'enable', 'on');
        
        set(handles.edit_freq_sweep, 'enable', 'off');
        set(handles.edit_layer_number, 'enable', 'off');
        set(handles.popupmenu_rt, 'enable', 'off');
end

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_sweep contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_sweep


% --- Executes during object creation, after setting all properties.
function popupmenu_sweep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_sweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_solution_type.
function popupmenu_solution_type_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_solution_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_solution_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_solution_type


% --- Executes during object creation, after setting all properties.
function popupmenu_solution_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_solution_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_freq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_freq as text
%        str2double(get(hObject,'String')) returns contents of edit_freq as a double


% --- Executes during object creation, after setting all properties.
function edit_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_theta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_theta as text
%        str2double(get(hObject,'String')) returns contents of edit_theta as a double


% --- Executes during object creation, after setting all properties.
function edit_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_freq_sweep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_freq_sweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_freq_sweep as text
%        str2double(get(hObject,'String')) returns contents of edit_freq_sweep as a double


% --- Executes during object creation, after setting all properties.
function edit_freq_sweep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_freq_sweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_rt.
function popupmenu_rt_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_rt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_rt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_rt


% --- Executes during object creation, after setting all properties.
function popupmenu_rt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_rt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_layer_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_layer_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_layer_number as text
%        str2double(get(hObject,'String')) returns contents of edit_layer_number as a double


% --- Executes during object creation, after setting all properties.
function edit_layer_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_layer_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_layer_number2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_layer_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_layer_number as text
%        str2double(get(hObject,'String')) returns contents of edit_layer_number as a double


% --- Executes during object creation, after setting all properties.
function edit_layer_number2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_layer_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_rt.
function popupmenu_rt2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_rt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_rt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_rt


% --- Executes during object creation, after setting all properties.
function popupmenu_rt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_rt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_theta_sweep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_theta_sweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_theta_sweep as text
%        str2double(get(hObject,'String')) returns contents of edit_theta_sweep as a double


% --- Executes during object creation, after setting all properties.
function edit_theta_sweep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_theta_sweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_bc.
function popupmenu_bc_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_bc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_bc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_bc


% --- Executes during object creation, after setting all properties.
function popupmenu_bc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_bc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_calculate.
function pushbutton_calculate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Ez

freq = eval(get(handles.edit_freq, 'string'))*1e9;
inc_ang = eval(get(handles.edit_theta, 'string'));

eps_r = eval(get(handles.edit_eps_r, 'string'));
mu_r = eval(get(handles.edit_mu_r, 'string'));
thickness = eval(get(handles.edit_thickness, 'string'))*1e-3;

eps_0 = 8.854187e-12;
mu_0 = 4*pi*1e-7;

k0 = 2*pi*freq*sqrt(eps_0 * mu_0);

minC = str2double(get(handles.edit_minC, 'string'));
maxC = str2double(get(handles.edit_maxC, 'string'));

x = eval(get(handles.edit_x, 'string'))*1e-3;
y = eval(get(handles.edit_y, 'string'))*1e-3;

bc = get(handles.popupmenu_bc, 'value');
solution_type = get(handles.popupmenu_solution_type, 'value');

switch bc 
    case 1
        yy_end = 0;
    case 2
        yy_end = -1;
    case 3
        yy_end = 1;
    case 4
        yy_end = 0;
end

 hSelectedDirection = get(handles.uibuttongroup2, 'SelectedObject');
  selectedDirection = get(hSelectedDirection, 'Tag');
  
  
%Construct Layer
minY  = y(1);
maxY = y(end);

nX = length(x);
nY = length(y);

nnLayer = length(eps_r);

set(handles.axes2, 'visible', 'on');
axes(handles.axes2);

axis off;

for ii = 1:nnLayer
    
    txt_x_location = (sum(thickness(1:ii))+sum(thickness(1:(ii-1))))/2*1e3;
    txt_y_location = 0;
    
    txt = ['Layer', num2str(ii)];
    txt_eps = ['\epsilon_r = ', num2str(eps_r(ii))];
    txt_mu = ['\mu_r = ', num2str(mu_r(ii))];
    
    r = 0;
    g = 0;
    b = 0;
    
    text(txt_x_location, txt_y_location, txt, 'color', [r g b], 'fontsize', 15);
    
end

d = zeros(size(thickness));
kd = zeros(size(thickness));
%Propagation Constant of i'th Later
for ii = 1:length(d)
    
    d(ii) = sum(thickness(1:ii));
    
    if (eps_r(ii) < 0 && mu_r(ii) < 0)
        %Metamaterial (Double Negative Media)
        kd(ii) = -k0*sqrt(abs(eps_r(ii)*abs(mu_r(ii))));
    else
        kd(ii) = k0*sqrt(eps_r(ii)*mu_r(ii));
    end
    
end

theta = zeros(size(kd));


%Phase Matching Condition (Propagation Angle of 1'th Layer)
for ii = 1:length(kd)
    
    if ii == 1
        theta(ii) = inc_ang;
    else
        theta(ii) = asind(      kd(ii - 1)/kd(ii)     *    sind(theta(ii-1))    );
    end
    
end

%Coefficients Calculation
for ii = nnLayer:-1:2
    
    if ii == nnLayer
        yy(ii - 1) = recursive_model(yy_end, d(ii - 1), kd(ii-1)*cosd(theta(ii - 1)), kd(ii)*cosd(theta(ii)), mu_0*mu_r(ii-1), mu_0*mu_r(ii));
    else
        yy(ii - 1) = recursive_model(yy(ii), d(ii - 1), kd(ii-1)*cosd(theta(ii - 1)), kd(ii)*cosd(theta(ii)), mu_0*mu_r(ii-1), mu_0*mu_r(ii));
    end
    
end

%Calculation of Magnitudes of Fields
mag_plus = zeros(size(thickness));
mag_minus = zeros(size(thickness));

for ii = 1:nnLayer
    
    if ii == 1
        mag_plus(ii) = 1;
        mag_minus(ii) = mag_plus(ii)*yy(ii);
    else
        [mag_plus(ii), mag_minus(ii)] = forward_calculation(mag_plus(ii - 1), mag_minus(ii - 1), d(ii - 1), kd(ii-1)*cosd(theta(ii - 1)), kd(ii)*cosd(theta(ii)), mu_0*mu_r(ii-1), mu_0*mu_r(ii));
    end
    
end


Ez = plot_wave(thickness, d , kd, theta, mag_plus, mag_minus, minY,maxY, minC, maxC, nX, nY, nnLayer, 2*pi*freq,selectedDirection , hObject,handles);








% --- Executes on button press in pushbutton_animate.
function pushbutton_animate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_animate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Ez

thickness = eval(get(handles.edit_thickness, 'string'));

nFrame = 100;
duration = 10;

animate_field(Ez, thickness, nFrame, duration, handles);



function edit_minC_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minC as text
%        str2double(get(hObject,'String')) returns contents of edit_minC as a double


% --- Executes during object creation, after setting all properties.
function edit_minC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxC_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxC as text
%        str2double(get(hObject,'String')) returns contents of edit_maxC as a double


% --- Executes during object creation, after setting all properties.
function edit_maxC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_x as text
%        str2double(get(hObject,'String')) returns contents of edit_x as a double

x = eval(get(handles.edit_x, 'string'));
thickness = eval(get(handles.edit_thickness, 'string'));

max_X_thickness = sum(thickness);

if (x(end) > max_X_thickness)
    
    extra_x = x(end) - max_X_thickness;
    thickness(end) = thickness(end) + extra_x;
    set(handles.edit_thickness, 'string', mat2str(thickness));
    
end




% --- Executes during object creation, after setting all properties.
function edit_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_y as text
%        str2double(get(hObject,'String')) returns contents of edit_y as a double


% --- Executes during object creation, after setting all properties.
function edit_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_prop_vec.
function popupmenu_prop_vec_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_prop_vec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_prop_vec contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_prop_vec


% --- Executes during object creation, after setting all properties.
function popupmenu_prop_vec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_prop_vec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_thickness_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ mu_r, eps_r] = eval_mu_r_eps_r_from_thickness(handles);

set(handles.edit_eps_r, 'string', mat2str(eps_r));
set(handles.edit_mu_r, 'string', mat2str(mu_r));

thickness = eval(get(handles.edit_thickness, 'string'));
x = eval(get(handles.edit_x, 'string'));

delta = x(2) - x(1);
delta_str = num2str(delta);

maxX = sum(thickness);

maxX_str = num2str(maxX);

minX = x(1);

minX_str = num2str(minX);

x_str = strcat(minX_str, ':', delta_str, ':', maxX_str);

set(handles.edit_x, 'string', x_str);



% Hints: get(hObject,'String') returns contents of edit_thickness as text
%        str2double(get(hObject,'String')) returns contents of edit_thickness as a double


% --- Executes during object creation, after setting all properties.
function edit_thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_construct_model.
function pushbutton_construct_model_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_construct_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

eps_r = eval(get(handles.edit_eps_r, 'string'));
mu_r = eval(get(handles.edit_mu_r, 'string'));
thickness = eval(get(handles.edit_thickness, 'string'))*1e-3;

axes(handles.axes2);
cla;
set(handles.axes2, 'visible', 'off');


nnLayer = length(eps_r);

x = eval(get(handles.edit_x, 'string'))*1e-3;
y = eval(get(handles.edit_y, 'string'))*1e-3;

minX = x(1);
maxX = x(end);

minY = y(1);
maxY = y(end);

axes(handles.axes1);
cla;
axis([minX, maxX, minY, maxY]*1e3);
colorbar off;

set(handles.axes1, 'YDir', 'Normal');

for ii = 1:nnLayer
    
    x_location = sum(thickness(1:ii));
    txt_y_location = 0;
    txt_x_location = (sum(thickness(1:ii))+sum(thickness(1:(ii-1))))/2*1e3;
    
    txt = ['Layer', num2str(ii)];
    txt_eps = ['\epsilon_r = ', num2str(eps_r(ii))];
    txt_mu = ['\mu_r = ', num2str(mu_r(ii))];
    r = 0;
    g = 0;
    b = 0;
    
    if (x_location > minX && x_location < maxX)
        
        line([x_location x_location]*1e3, [minY, maxY]*1e3, 'color', 'black', 'linestyle', '--', 'linewidth', 2);

    end
    
    text(txt_x_location, txt_y_location, txt, 'color', [r g b], 'fontsize', 10);
    text(txt_x_location, txt_y_location-20,txt_eps, 'fontsize', 10, 'color', [r g b], 'fontweight', 'bold' );
    text(txt_x_location, txt_y_location-40,txt_mu, 'fontsize', 10, 'color', [r g b], 'fontweight', 'bold' );
    
end




% --- Executes on button press in pushbutton_clear.
function pushbutton_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);

cla;
colorbar off;

axes(handles.axes2);

cla;
colorbar off;



function edit_eps_r_Callback(hObject, eventdata, handles)
% hObject    handle to edit_eps_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[thickness, mu_r] = eval_thickness_mu_r_from_eps_r(handles);

set(handles.edit_thickness, 'string', mat2str(thickness));
set(handles.edit_mu_r, 'string', mat2str(mu_r));

% Hints: get(hObject,'String') returns contents of edit_eps_r as text
%        str2double(get(hObject,'String')) returns contents of edit_eps_r as a double


% --- Executes during object creation, after setting all properties.
function edit_eps_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_eps_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mu_r_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mu_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mu_r as text
%        str2double(get(hObject,'String')) returns contents of edit_mu_r as a double

[eps_r, thickness] = eval_eps_r_thickness_from_mu_r(handles)

set(handles.edit_thickness, 'string', mat2str(thickness));
set(handles.edit_eps_r, 'string', mat2str(eps_r));


% --- Executes during object creation, after setting all properties.
function edit_mu_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mu_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alphaX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alphaX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alphaX as text
%        str2double(get(hObject,'String')) returns contents of edit_alphaX as a double


% --- Executes during object creation, after setting all properties.
function edit_alphaX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alphaX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alphaY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alphaY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alphaY as text
%        str2double(get(hObject,'String')) returns contents of edit_alphaY as a double


% --- Executes during object creation, after setting all properties.
function edit_alphaY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alphaY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_betaX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_betaX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_betaX as text
%        str2double(get(hObject,'String')) returns contents of edit_betaX as a double


% --- Executes during object creation, after setting all properties.
function edit_betaX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_betaX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_betaY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_betaY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_betaY as text
%        str2double(get(hObject,'String')) returns contents of edit_betaY as a double


% --- Executes during object creation, after setting all properties.
function edit_betaY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_betaY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function varargout = NetworkViewer(varargin)
% NETWORKVIEWER MATLAB code for NetworkViewer.fig
%      NETWORKVIEWER, by itself, creates a new NETWORKVIEWER or raises the existing
%      singleton*.
%
%      H = NETWORKVIEWER returns the handle to a new NETWORKVIEWER or the handle to
%      the existing singleton*.
%
%      NETWORKVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NETWORKVIEWER.M with the given input arguments.
%
%      NETWORKVIEWER('Property','Value',...) creates a new NETWORKVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NetworkViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NetworkViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NetworkViewer

% Last Modified by GUIDE v2.5 05-Dec-2017 10:11:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NetworkViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @NetworkViewer_OutputFcn, ...
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


% --- Executes just before NetworkViewer is made visible.
function NetworkViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NetworkViewer (see VARARGIN)

% Choose default command line output for NetworkViewer
handles.output = hObject;


%*****************************************************************
%Change Read In data to be provided by NetworkViewer function call
%*****************************************************************

%Handle user argument for network type
if nargin > 4
    %Read in connectivity data
    switch(varargin{1})
        case 'coh'
            handles.data = evalin('base','ECoG_conn.coh'); 
        case 'wpli'
            handles.data = evalin('base','ECoG_conn.WPLI'); 
        otherwise
            error('Please provide a connectivity metric. Options are ''coh'' or ''wpli''')
    end
    switch(varargin{2})
        case 'unweighted'
            handles.NetworkType = 'Binarized';
        case 'weighted'
            handles.NetworkType = 'Thresholded';
        case 'u'
            handles.NetworkType = 'Binarized';
        case 'w'
            handles.NetworkType = 'Thresholded';
        otherwise
            error('Please provide a network type. Options are ''weighted'' or ''unweighted''. Defaulting to ''unweighted''')  
    end
else
    error('Please provide the following arguments: [Connectivity Type] [Network Type]')    
end

%Read in patient data
handles.patientdata = evalin('base','batchParams');


%Initialize popup boxes
set(handles.PatientMenu,'String',fieldnames(handles.data));
Patients = cellstr(get(handles.PatientMenu,'String'));
currentPatient = Patients{get(handles.PatientMenu,'Value')};
set(handles.ConditionMenu,'String',fieldnames(handles.data.(currentPatient)));
Conditions = cellstr(get(handles.ConditionMenu,'String'));
currentCondition = Conditions{get(handles.ConditionMenu,'Value')};
set(handles.BandMenu,'String',fieldnames(handles.data.(currentPatient).(currentCondition)));
Bands = cellstr(get(handles.BandMenu,'String'));
currentBand = Bands{get(handles.BandMenu,'Value')};


%Update string displaying network type
if strcmp(handles.NetworkType, 'Thresholded')
    set(handles.NetTypeString, 'String', 'Network Type:   Weighted');
elseif strcmp(handles.NetworkType, 'Binarized')
    set(handles.NetTypeString, 'String', 'Network Type:   Unweighted');
end

%Initialize first percolation threshold string
set(handles.PercThrString, 'String', ['Percolation Threshold:   ', num2str(handles.data.(currentPatient).(currentCondition).(currentBand).PercThr)]);

%Initialize Flags
handles.CoordFlag = 0;
handles.NodeLabelFlag = 0;

%Initialize Data Table
netType = handles.NetworkType;
modularity = handles.data.(currentPatient).(currentCondition).(currentBand).NetworkStats.(netType).Modularity;
transit = handles.data.(currentPatient).(currentCondition).(currentBand).NetworkStats.(netType).Transitivity;
charpathlen = handles.data.(currentPatient).(currentCondition).(currentBand).NetworkStats.(netType).CharPathLen;
swp = handles.data.(currentPatient).(currentCondition).(currentBand).NetworkStats.(netType).SWP;
set(handles.DataTable,'Data',[modularity,transit(:,2),charpathlen(:,2),swp(:,2)])
set(handles.DataTable,'ColumnName',{'Thr','Modularity','Trans','CharPath','Small-World'})


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NetworkViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NetworkViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in PatientMenu.
function PatientMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PatientMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PatientMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PatientMenu
Patients = cellstr(get(hObject,'String'));
currentPatient = Patients{get(hObject,'Value')};

%Update ConditionMenu with conditons available for selected patient
set(handles.ConditionMenu,'String',fieldnames(handles.data.(currentPatient)));

%Update BandMenu with bands avaiable for selected conditon
Conditions = cellstr(get(handles.ConditionMenu,'String'));
currentCondition = Conditions{1};
set(handles.BandMenu,'String',fieldnames(handles.data.(currentPatient).(currentCondition)));



% --- Executes during object creation, after setting all properties.
function PatientMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PatientMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in ConditionMenu.
function ConditionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ConditionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ConditionMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ConditionMenu


% --- Executes during object creation, after setting all properties.
function ConditionMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ConditionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BandMenu.
function BandMenu_Callback(hObject, eventdata, handles)
% hObject    handle to BandMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BandMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BandMenu


% --- Executes during object creation, after setting all properties.
function BandMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BandMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ModularityPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModularityPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ModularityPlot


% --- Executes during object creation, after setting all properties.
function TransPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TransPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate TransPlot



% --- Executes on button press in UpdateButton.
function UpdateButton_Callback(hObject, eventdata, handles)
% hObject    handle to UpdateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Get the current selected state for plotting
Patients = cellstr(get(handles.PatientMenu,'String'));
Conditions = cellstr(get(handles.ConditionMenu,'String'));
Bands = cellstr(get(handles.BandMenu,'String'));
currentPatient = Patients{get(handles.PatientMenu,'Value')};
currentCondition = Conditions{get(handles.ConditionMenu,'Value')};
currentBand = Bands{get(handles.BandMenu,'Value')};
netType = handles.NetworkType;

%Load graph metrics for selected data
perc = handles.data.(currentPatient).(currentCondition).(currentBand).PercThr;
modularity = handles.data.(currentPatient).(currentCondition).(currentBand).NetworkStats.(netType).Modularity;
transit = handles.data.(currentPatient).(currentCondition).(currentBand).NetworkStats.(netType).Transitivity;
charpathlen = handles.data.(currentPatient).(currentCondition).(currentBand).NetworkStats.(netType).CharPathLen;
swp = handles.data.(currentPatient).(currentCondition).(currentBand).NetworkStats.(netType).SWP;

%Update percolation threshold string
set(handles.PercThrString, 'String', ['Percolation Threshold:   ', num2str(perc)]);

%Plot Modularity
axes(handles.ModularityPlot);
if(isfield(handles,'modularity_l') && any(strcmp(handles.modularity_l.String, [currentPatient,' ',currentCondition,' ', currentBand])))
else
    modularity_p = plot(modularity(:,1),modularity(:,2),'-o','DisplayName',[currentPatient,' ',currentCondition,' ', currentBand]);
    hold on;
    modularity_perc = scatter(modularity(find(modularity(:,1) == perc),1),modularity(find(modularity(:,1) == perc),2),'*r','DisplayName','Percolation Threshold');
    handles.modularity_l = legend('-DynamicLegend','Location','northwest');
    set(gca, 'Xdir', 'reverse')
    xlabel('Percent of Edges Retained');
    ylabel('Modularity');
    guidata(hObject, handles);   
end


%Plot Transitivity
axes(handles.TransPlot);
if(isfield(handles,'transit_l') && any(strcmp(handles.transit_l.String, [currentPatient,' ',currentCondition,' ', currentBand])))
else
    transit_p = plot(transit(:,1),transit(:,2),'-o','DisplayName',[currentPatient,' ',currentCondition,' ', currentBand]);
    hold on;
    transit_perc = scatter(transit(find(transit(:,1) == perc),1),transit(find(transit(:,1) == perc),2),'*r','DisplayName','Percolation Threshold');
    handles.transit_l = legend('-DynamicLegend','Location','northwest');
    set(gca, 'Xdir', 'reverse')
    xlabel('Percent of Edges Retained');
    ylabel('Transitivity');
    guidata(hObject, handles); 
end



%Plot Characteristic Path Length
axes(handles.CharPathPlot);
if(isfield(handles,'charpathlen_l') && any(strcmp(handles.charpathlen_l.String, [currentPatient,' ',currentCondition,' ', currentBand])))
else
    charpathlen_p = plot(charpathlen(:,1),charpathlen(:,2),'-o','DisplayName',[currentPatient,' ',currentCondition,' ', currentBand]);
    hold on;
    charpathlen_perc = scatter(charpathlen(find(charpathlen(:,1) == perc),1),charpathlen(find(charpathlen(:,1) == perc),2),'*r','DisplayName','Percolation Threshold');
    handles.charpathlen_l = legend('-DynamicLegend','Location','northwest');
    set(gca, 'Xdir', 'reverse')
    xlabel('Percent of Edges Retained');
    ylabel('Characteristic Path Length');
    guidata(hObject, handles); 
end


%Plot Small-World Propensity
axes(handles.SWPlot);
if(isfield(handles,'swp_l') && any(strcmp(handles.swp_l.String, [currentPatient,' ',currentCondition,' ', currentBand])))
else
    swp_p = plot(swp(:,1),swp(:,2),'-o','DisplayName',[currentPatient,' ',currentCondition,' ', currentBand]);
    hold on;
    swp_perc = scatter(swp(find(swp(:,1) == perc),1),swp(find(swp(:,1) == perc),2),'*r','DisplayName','Percolation Threshold');
    handles.swp_l = legend('-DynamicLegend','Location','northwest');
    set(gca, 'Xdir', 'reverse')
    xlabel('Percent of Edges Retained');
    ylabel('Small-World Propensity');
    guidata(hObject, handles);
end

%Update Data Table
set(handles.DataTable,'Data',[modularity,transit(:,2),charpathlen(:,2),swp(:,2)])


% --- Executes on mouse press over axes background.
function ModularityPlot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ModularityPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function DataTable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected cell(s) is changed in DataTable.
function DataTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to DataTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

selectedRow = eventdata.Indices(1);
data=get(hObject,'Data');
data=data(selectedRow,1);


UpdateNetwork(data, handles);

function UpdateNetwork(data, handles)
%Get data on current state
Patients = cellstr(get(handles.PatientMenu,'String'));
Conditions = cellstr(get(handles.ConditionMenu,'String'));
Bands = cellstr(get(handles.BandMenu,'String'));
currentPatient = Patients{get(handles.PatientMenu,'Value')};
currentCondition = Conditions{get(handles.ConditionMenu,'Value')};
currentBand = Bands{get(handles.BandMenu,'Value')};
netType = handles.NetworkType;


%Get electrode data for current patient
channels = handles.patientdata.(currentPatient).(currentCondition).ECoGchannels;
xCoord = [channels.xCoord];
yCoord = [channels.yCoord];
zCoord = [channels.zCoord];
nodelabs = {channels.ROI};

%Make network plot for threshold of selected row
Networks = fieldnames(handles.data.(currentPatient).(currentCondition).(currentBand).Networks.(netType));
Networks = [Networks, strrep(Networks,'thr','')];
netIndex = find(strcmp(Networks(:,2),num2str(data)));
if(isempty(netIndex))
    netIndex = find(strcmp(Networks(:,2),'prc'));
end
currentNet = Networks{netIndex};
G = graph(handles.data.(currentPatient).(currentCondition).(currentBand).Networks.(netType).(currentNet));
axes(handles.NetPlot);

if(handles.CoordFlag == 1)
    if(handles.NodeLabelFlag == 1)
        plot(G,'XData',xCoord,'YData',yCoord, 'NodeLabel', nodelabs)
    elseif(handles.NodeLabelFlag == 0)
        plot(G,'XData',xCoord,'YData',yCoord)
    end
elseif(handles.CoordFlag == 0)
    if(handles.NodeLabelFlag == 1)
        plot(G, 'NodeLabel', nodelabs)
    elseif(handles.NodeLabelFlag == 0)
        plot(G)
    end
end




% --- Executes on button press in CoordCheck.
function CoordCheck_Callback(hObject, eventdata, handles)
% hObject    handle to CoordCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CoordCheck
if (get(hObject,'Value') == get(hObject,'Max'))
	handles.CoordFlag = 1;
else
	handles.CoordFlag = 0;
end
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in LabelCheck.
function LabelCheck_Callback(hObject, eventdata, handles)
% hObject    handle to LabelCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LabelCheck
if (get(hObject,'Value') == get(hObject,'Max'))
	handles.NodeLabelFlag = 1;
else
	handles.NodeLabelFlag = 0;
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function NetPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NetPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate NetPlot


% --- Executes on button press in ClearButton.
function ClearButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ModularityPlot);
cla
handles = rmfield(handles, 'modularity_l');
guidata(hObject, handles);
legend('hide');

axes(handles.TransPlot);
cla
handles = rmfield(handles, 'transit_l');
guidata(hObject, handles);
legend('hide');

axes(handles.CharPathPlot);
cla
handles = rmfield(handles, 'charpathlen_l');
guidata(hObject, handles);
legend('hide');

axes(handles.SWPlot);
cla
handles = rmfield(handles, 'swp_l');
guidata(hObject, handles);
legend('hide');


% --- Executes during object creation, after setting all properties.
function SWPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SWPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate SWPlot

function varargout = flag_diag(varargin)
% FLAG_DIAG MATLAB code for flag_diag.fig
%      FLAG_DIAG, by itself, creates a new FLAG_DIAG or raises the existing
%      singleton*.
%
%      H = FLAG_DIAG returns the handle to a new FLAG_DIAG or the handle to
%      the existing singleton*.
%
%      FLAG_DIAG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLAG_DIAG.M with the given input arguments.
%
%      FLAG_DIAG('Property','Value',...) creates a new FLAG_DIAG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before flag_diag_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to flag_diag_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help flag_diag

% Last Modified by GUIDE v2.5 15-Feb-2018 15:54:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @flag_diag_OpeningFcn, ...
    'gui_OutputFcn',  @flag_diag_OutputFcn, ...
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


% --- Executes just before flag_diag is made visible.
function flag_diag_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to flag_diag (see VARARGIN)
n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);
%keyboard
if isfield(s,'Config')==1
    handles.floatname=s.Config.floatname;
    handles.dacname=s.Config.dacname;
    handles.DIR_FTP=s.Config.DIR_FTP;
    handles.FILE_TOPO=s.Config.FILE_TOPO;
	handles.VPN=s.Config.VPN;
	handles.FLAG=s.Config.FLAG;
else
    disp('Warning:********* enter CONFIG as un imput argument: flag_diag(''Config'',CONFIG)')
    handles.floatname='';
    handles.dacname='';
    handles.DIR_FTP='';
	handles.VPN='';
	handles.FLAG;
end
%keyboard
set(hObject,'toolbar','Figure');
% data for the plot
set(handles.figure1,'Color',[1 1 1]);
[handles.file_list]=select_float_files_on_ftp(handles.floatname,handles.dacname,handles.DIR_FTP,'C');
[handles.file_listB]=select_float_files_on_ftp(handles.floatname,handles.dacname,handles.DIR_FTP,'B');

% Choose default command line output for flag_diag
handles.output = hObject;
set(handles.profile,'String','');


Topo=read_netcdf_allthefile(handles.FILE_TOPO);

handles.topo=Topo;
handles.ax1=subplot(4,18,[1:5,19:23]);
hold on
box on
handles.ax2=subplot(1,12,[5:9]);
box on
handles.ax3=subplot(4,18,[37:38,55:56]);
hold on
box on
handles.ax4=subplot(4,18,[40:41,58:59]);
hold on
box on
% Update handles structure
guidata(hObject, handles);


%keyboard
% UIWAIT makes flag_diag wait for user response (see UIRESUME)
 %uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = flag_diag_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;



function new_flag=new_flag_Callback(hObject, eventdata, handles)
% hObject    handle to new_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of new_flag as text
%        str2double(get(hObject,'String')) returns contents of new_flag as a double
parameters=get(handles.parameter,'String');
choice=get(handles.parameter,'Value');
if strcmp(parameters{choice},'Temperature')||strcmp(parameters{choice},'Salinity')||strcmp(parameters{choice},'ThetaS')||strcmp(parameters{choice},'Oxygene')
new_flag=str2double(get(hObject,'String'));
BackgroundColor_vec=[0 1 0; 1 1 0; 1 0 1; 1 0 0];
set(hObject,'BackgroundColor',BackgroundColor_vec(new_flag,:))
else
    helpdlg(['Flags cannot be changed for ' parameters{choice}])
end
    
% --- Executes during object creation, after setting all properties.
function new_flag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to new_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selection.
function selection_Callback(hObject, eventdata, handles)
% hObject    handle to selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%keyboard
parameters=get(handles.selection,'String');
choice=get(handles.selection,'Value');

new_flag=get(handles.new_flag,'String');
colormarker={'g','y','m','r'};
xx=[];
yy=[];

%zoom out
zoom off
if handles.VPN==0
	if ~isempty(strfind(parameters{choice},'Rectangle'))
		[h,x,y]=selectdata('selectionmode','Rect','Pointer','cross')
	end
	if ~isempty(strfind(parameters{choice},'Single Point'))
		[h,x,y]=selectdata('selectionmode','Closest','Pointer','cross')

	end

	for l=1:length(x)
		xx=[xx x{l}']; yy=[yy y{l}'];
	end
	min_x=min(xx);
	min_y=min(yy);
	max_x=max(xx);
	max_y=max(yy);
end
if handles.VPN==1
	if ~isempty(strfind(parameters{choice},'Rectangle'))
		
		[xx1,yy1]=ginput(1);plot(xx1,yy1,'+')
		[xx2,yy2]=ginput(1);plot([xx1 xx2],[yy1 yy2],'-+')
		[xx3,yy3]=ginput(1);plot([xx1 xx2 xx3],[yy1 yy2 yy3],'-+')
		[xx4,yy4]=ginput(1);plot([xx1 xx2 xx3 xx4 xx1],[yy1 yy2 yy3 yy4 yy1],'-+')
		xx=[xx1 xx2 xx3 xx4];
		yy=[yy1 yy2 yy3 yy4];
		min_x=min(xx);
		min_y=min(yy);
		max_x=max(xx);
		max_y=max(yy);
	end

	if ~isempty(strfind(parameters{choice},'Single Point'))
		[h,x,y]=selectdata('selectionmode','Closest','Pointer','cross');
		for l=1:length(x)
			xx=[xx x{l}']; yy=[yy y{l}'];
		end
		min_x=min(xx);
		min_y=min(yy);
		max_x=max(xx);
		max_y=max(yy);
	end
end
parameters=get(handles.parameter,'String');
choice=get(handles.parameter,'Value');
prof_number=get(handles.profile,'String');

if isempty(prof_number)
    set(handles.profile,'String','1');
    profile_Callback(hObject, eventdata, handles)
else
    [d,r]=strtok(prof_number,':');d=str2num(d);
    if isempty(r)==0
        r=strrep(r,'end',num2str(size(handles.F.cycle_number.data,1)));
        r=str2num(strrep(r,':',''));
    else
        r=d;
    end
    %keyboard
    for k=1:length(handles.thelist)
%     for k=1:r-d+1
    if ~isempty(min_x)|~isempty(max_x)|~isempty(min_y)|~isempty(max_y)
        if strcmp(parameters{choice},'Temperature')
            isub=find(handles.F.temp.data(k,:)<=max_x&handles.F.temp.data(k,:)>=min_x&handles.F.pres.data(k,:)<=max_y&handles.F.pres.data(k,:)>=min_y);
            hold on
            %keyboard
            plot(handles.F.temp.data(k,isub),handles.F.pres.data(k,isub),[colormarker{str2num(new_flag)} '+'])
            %handles.F.temp_qc_new=handles.F.temp_qc;
            handles.F.temp_qc.data(k,isub)=str2num(new_flag);
        elseif strcmp(parameters{choice},'Salinity')
            isub=find(handles.F.psal.data(k,:)<=max_x&handles.F.psal.data(k,:)>=min_x&handles.F.pres.data(k,:)<=max_y&handles.F.pres.data(k,:)>=min_y);
            hold on
            
            plot(handles.F.psal.data(k,isub),handles.F.pres.data(k,isub),[colormarker{str2num(new_flag)} '+'])
            % handles.F.psal_qc_new=handles.F.psal_qc;
            handles.F.psal_qc.data(k,isub)=str2num(new_flag);
            
        elseif strcmp(parameters{choice},'ThetaS')
            ptemp = sw_ptmp(handles.F.psal.data,handles.F.temp.data,handles.F.pres.data,0);
            isub=find(handles.F.psal.data(k,:)<=max_x&handles.F.psal.data(k,:)>=min_x&ptemp(k,:)<=max_y&ptemp(k,:)>=min_y);
            hold on
            
            plot(handles.F.psal.data(k,isub),ptemp(k,isub),[colormarker{str2num(new_flag)} '+'])
            % handles.F.psal_qc_new=handles.F.psal_qc;
            handles.F.psal_qc.data(k,isub)=str2num(new_flag);    
		elseif strcmp(parameters{choice},'Oxygene')
            isub=find(handles.F.doxy.data(k,:)<=max_x&handles.F.doxy.data(k,:)>=min_x&handles.F.pres.data(k,:)<=max_y&handles.F.pres.data(k,:)>=min_y);
            hold on
            
            plot(handles.F.doxy.data(k,isub),handles.F.pres.data(k,isub),[colormarker{str2num(new_flag)} '+'])
            % handles.F.psal_qc_new=handles.F.psal_qc;
            handles.F.doxy_qc.data(k,isub)=str2num(new_flag);
            
%         elseif strcmp(parameters{choice},'Temperature(Adjusted)')
%             isub=find(handles.F.temp_adjusted.data(k,:)<=max_x&handles.F.temp_adjusted.data(k,:)>=min_x&handles.F.pres_adjusted.data(k,:)<=max_y&handles.F.pres_adjusted.data(k,:)>=min_y);
%             hold on
%             %keyboard
%             plot(handles.F.temp_adjusted.data(k,isub),handles.F.pres_adjusted.data(k,isub),[colormarker{str2num(new_flag)} '+'])
%             %handles.F.temp_qc_new=handles.F.temp_qc;
%             handles.F.temp_adjusted_qc.data(k,isub)=str2num(new_flag);
%         elseif strcmp(parameters{choice},'Salinity(Adjusted)')
%             isub=find(handles.F.psal.data(k,:)<=max_x&handles.F.psal.data(k,:)>=min_x&handles.F.pres.data(k,:)<=max_y&handles.F.pres.data(k,:)>=min_y);
%             hold on
%             
%             plot(handles.F.psal_adjusted.data(k,isub),handles.F.pres_adjusted.data(k,isub),[colormarker{str2num(new_flag)} '+'])
%             % handles.F.psal_qc_new=handles.F.psal_qc;
%             handles.F.psal_adjusted_qc.data(k,isub)=str2num(new_flag);
        end
    end
    end
end
set(handles.selection,'Value',1);
set(gca,'XLimMode','auto')
set(gca,'YLimMode','auto')
% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in parameter.
function parameter_Callback(hObject, eventdata, handles)
% hObject    handle to parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns parameter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parameter
parameters=get(handles.parameter,'String');
choice=get(handles.parameter,'Value');
%disp('je passe dans Parameter callback1')
if isfield(handles,'choice_previous')==1
 handles.choice_previous
 choice
	if handles.choice_previous~=6&&handles.choice_previous~=11
	   if choice==6||choice==11
	      handles.choice_previous=choice;
		  guidata(hObject,handles);
		  %disp('je passe dans Parameter callback1a')
		  profile_Callback(hObject, eventdata, handles)
		  return
		  
		end
	else
		if choice~=6&&choice~=11
		 handles.choice_previous=choice;
		 guidata(hObject,handles);
		 %disp('je passe dans Parameter callback1b')
		 profile_Callback(hObject, eventdata, handles)
		 return
		   
		end
	end
		
end
%disp('je passe dans Parameter callback2')
prof_number=get(handles.profile,'String');
zoom out
if isempty(prof_number)
    set(handles.profile,'String','1');
    profile_Callback(hObject, eventdata, handles)
	return
else
    axes(handles.ax2)
	cla
	axes(handles.ax1)
	if handles.FLAG==0
	   markerstr='none';
	else
	   markerstr='+';
	end
	markerstr
    	if isempty(handles.F.latitude.data)==0
	    cla
		[t,r]=strtok(handles.file_list,'_');
		r=strrep(r,'_','');
		r=strrep(r,'.nc','');
		r=strrep(r,'D','');
	   
		handles.F.min_cy=min(str2num(char(r)));
		handles.F.max_cy=max(str2num(char(r)));
		
		plot_traj_glob_gui(handles.F,handles.topo)
		if isempty(handles.F2)==0
		handles.F2.min_cy=min(str2num(char(r)));
		handles.F2.max_cy=max(str2num(char(r)));   
		plot_traj_glob_gui(handles.F2,handles.topo)
		end
		axes(handles.ax2)
		 
		%handles.F.psal.data
		parameters{choice}
		
		if strcmp(parameters{choice},'Temperature')
		    
			cla(gca)
			zoom out
			n_prof=size(handles.F.temp.data,1);
			
			if isempty(handles.F2)==0
			n_prof2=size(handles.F2.temp.data,1);
			plot_profile_with_flag(handles.F2,'temp','pres',[1:n_prof2],'LineWidth',0.5,'Marker','none');
			end
			thetitle=plot_profile_with_flag(handles.F,'temp','pres',[1:n_prof],'LineWidth',0.5,'Marker',markerstr);
			title(thetitle)
			axes(handles.ax3)
			cla(gca)
			plot_profile_with_flag(handles.F3,'temp','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax4)
			cla(gca)
			plot_profile_with_flag(handles.F3,'psal','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax2)
		elseif strcmp(parameters{choice},'Salinity')
			
			cla(gca)
			zoom out
			n_prof=size(handles.F.temp.data,1);        
			if isempty(handles.F2)==0
			n_prof2=size(handles.F2.temp.data,1);
			
			plot_profile_with_flag(handles.F2,'psal','pres',[1:n_prof2],'LineWidth',0.5,'Marker','none');
			end
			%disp('test')
			%markerstr
			thetitle= plot_profile_with_flag(handles.F,'psal','pres',[1:n_prof],'LineWidth',0.5,'Marker',markerstr);
			%thetitle= plot_profile_with_flag(handles.F,'psal','pres',[1:n_prof],'LineWidth',0.5,'Marker','none');
			title(thetitle)
			axes(handles.ax3)
			cla(gca)
			plot_profile_with_flag(handles.F3,'temp','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax4)
			cla(gca)
			plot_profile_with_flag(handles.F3,'psal','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax2)
		elseif strcmp(parameters{choice},'Oxygene')
			
			cla(gca)
			zoom out
			n_prof=size(handles.F.pres.data,1);        
			if isempty(handles.F2)==0
			n_prof2=size(handles.F2.pres.data,1);
			plot_profile_with_flag(handles.F2,'doxy','pres',[1:n_prof2],'LineWidth',0.5,'Marker','none');
			end
			thetitle= plot_profile_with_flag(handles.F,'doxy','pres',[1:n_prof],'LineWidth',0.5,'Marker',markerstr);
			%thetitle= plot_profile_with_flag(handles.F,'psal','pres',[1:n_prof],'LineWidth',0.5,'Marker','none');
			title(thetitle)
			axes(handles.ax3)
			cla(gca)
			plot_profile_with_flag(handles.F3,'temp','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax4)
			cla(gca)
			plot_profile_with_flag(handles.F3,'psal','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax2)
		 elseif strcmp(parameters{choice},'Ptemp')
		   % keyboard
			cla(gca)
			zoom out
			n_prof=size(handles.F.temp.data,1);   
			if isempty(handles.F2)==0
			n_prof2=size(handles.F2.temp.data,1);
			handles.F2.ptemp.data=sw_ptmp(handles.F2.psal.data,handles.F2.temp.data,handles.F2.pres.data,0);
			handles.F2.ptemp_qc.data=handles.F2.psal_qc.data;
			plot_profile_with_flag(handles.F2,'ptemp','pres',[1:n_prof2],'LineWidth',0.5,'Marker','none');
			end
			handles.F.ptemp.data=sw_ptmp(handles.F.psal.data,handles.F.temp.data,handles.F.pres.data,0);
			handles.F.ptemp_qc.data=handles.F.psal_qc.data;
			thetitle= plot_profile_with_flag(handles.F,'ptemp','pres',[1:n_prof],'LineWidth',0.5,'Marker',markerstr);
			%thetitle= plot_profile_with_flag(handles.F,'ptemp','pres',[1:n_prof],'LineWidth',0.5,'Marker','none');
			title(thetitle)    
		    axes(handles.ax3)
			cla(gca)
			plot_profile_with_flag(handles.F3,'temp','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax4)
			cla(gca)
			plot_profile_with_flag(handles.F3,'psal','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax2)
		elseif strcmp(parameters{choice},'ThetaS')
		   % keyboard
			cla(gca)
			zoom out
			n_prof=size(handles.F.temp.data,1);      
			if isempty(handles.F2)==0
			n_prof2=size(handles.F2.temp.data,1);
			handles.F2.ptemp.data=sw_ptmp(handles.F2.psal.data,handles.F2.temp.data,handles.F2.pres.data,0);
			plot_profile_with_flag(handles.F2,'psal','ptemp',[1:n_prof2],'LineWidth',0.5,'Marker','none');
			end
			handles.F.ptemp=handles.F.temp;
			handles.F.ptemp.data=sw_ptmp(handles.F.psal.data,handles.F.temp.data,handles.F.pres.data,0);
			thetitle= plot_profile_with_flag(handles.F,'psal','ptemp',[1:n_prof],'LineWidth',0.5,'Marker',markerstr);
			title(thetitle)
			axes(handles.ax3)
			cla(gca)
			plot_profile_with_flag(handles.F3,'temp','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax4)
			cla(gca)
			plot_profile_with_flag(handles.F3,'psal','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax2)
		elseif strcmp(parameters{choice},'Temperature(Adjusted)')
			cla(gca)
			zoom out
			%keyboard
			n_prof=size(handles.F.temp.data,1);
			
			if isempty(handles.F2)==0
			n_prof2=size(handles.F2.temp.data,1);
			plot_profile_with_flag(handles.F2,'temp_adjusted','pres_adjusted',[1:n_prof2],'LineWidth',0.5,'Marker','none');
			end
			thetitle=plot_profile_with_flag(handles.F,'temp_adjusted','pres_adjusted',[1:n_prof],'LineWidth',0.5,'Marker',markerstr);
			title(thetitle)
			axes(handles.ax3)
			cla(gca)
			plot_profile_with_flag(handles.F3,'temp','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax4)
			cla(gca)
			plot_profile_with_flag(handles.F3,'psal','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax2)
		elseif strcmp(parameters{choice},'Ptemp(Adjusted)')
		   % keyboard
			cla(gca)
			zoom out
			n_prof=size(handles.F.temp.data,1);     
			if isempty(handles.F2)==0
			n_prof2=size(handles.F2.temp.data,1);
			handles.F2.ptemp_adjusted.data=sw_ptmp(handles.F2.psal_adjusted.data,handles.F2.temp_adjusted.data,handles.F2.pres_adjusted.data,0);
			handles.F2.ptemp_adjusted_qc.data=handles.F2.psal_adjusted_qc.data;
			plot_profile_with_flag(handles.F2,'ptemp_adjusted','pres_adjusted',[1:n_prof2],'LineWidth',0.5,'Marker','none');
			end
			handles.F.ptemp_adjusted.data=sw_ptmp(handles.F.psal_adjusted.data,handles.F.temp_adjusted.data,handles.F.pres_adjusted.data,0);
			handles.F.ptemp_adjusted_qc.data=handles.F.psal_adjusted_qc.data;
			thetitle= plot_profile_with_flag(handles.F,'ptemp_adjusted','pres_adjusted',[1:n_prof],'LineWidth',0.5,'Marker',markerstr);
			title(thetitle)   
			axes(handles.ax3)
			cla(gca)
			plot_profile_with_flag(handles.F3,'temp','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax4)
			cla(gca)
			plot_profile_with_flag(handles.F3,'psal','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax2)
		elseif strcmp(parameters{choice},'Salinity(Adjusted)')
			%keyboard
			cla(gca)
			zoom out
			n_prof=size(handles.F.temp.data,1);        
			if isempty(handles.F2)==0
			n_prof2=size(handles.F2.temp.data,1);
			plot_profile_with_flag(handles.F2,'psal_adjusted','pres_adjusted',[1:n_prof2],'LineWidth',0.5,'Marker','none');
			end
			thetitle= plot_profile_with_flag(handles.F,'psal_adjusted','pres_adjusted',[1:n_prof],'LineWidth',0.5,'Marker',markerstr);
			title(thetitle)
		   	axes(handles.ax3)
			cla(gca)
			plot_profile_with_flag(handles.F3,'temp','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax4)
			cla(gca)
			plot_profile_with_flag(handles.F3,'psal','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax2)
		 elseif strcmp(parameters{choice},'Oxygene(Adjusted)')
			
			cla(gca)
			zoom out
			n_prof=size(handles.F.pres.data,1);        
			if isempty(handles.F2)==0
			n_prof2=size(handles.F2.pres.data,1);
			plot_profile_with_flag(handles.F2,'doxy_adjusted','pres_adjusted',[1:n_prof2],'LineWidth',0.5,'Marker','none');
			end
			thetitle= plot_profile_with_flag(handles.F,'doxy_adjusted','pres_adjusted',[1:n_prof],'LineWidth',0.5,'Marker',markerstr);
			title(thetitle)
			axes(handles.ax3)
			cla(gca)
			plot_profile_with_flag(handles.F3,'temp','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax4)
			cla(gca)
			plot_profile_with_flag(handles.F3,'psal','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax2)
		elseif strcmp(parameters{choice},'ThetaS(Adjusted)')
		   % keyboard
			cla(gca)
			zoom out
			n_prof=size(handles.F.temp.data,1);  
			if isempty(handles.F2)==0
			n_prof2=size(handles.F2.temp.data,1);
			handles.F2.ptemp_adjusted.data=sw_ptmp(handles.F2.psal_adjusted.data,handles.F2.temp_adjusted.data,handles.F2.pres_adjusted.data,0);
			plot_profile_with_flag(handles.F2,'psal_adjusted','ptemp_adjusted',[1:n_prof2],'LineWidth',0.5,'Marker','none');
			end
			handles.F.ptemp_adjusted.data=sw_ptmp(handles.F.psal_adjusted.data,handles.F.temp_adjusted.data,handles.F.pres_adjusted.data,0);
			thetitle= plot_profile_with_flag(handles.F,'psal_adjusted','ptemp_adjusted',[1:n_prof],'LineWidth',0.5,'Marker',markerstr);
			title(thetitle) 
			axes(handles.ax3)
			cla(gca)
			plot_profile_with_flag(handles.F3,'temp','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax4)
			cla(gca)
			plot_profile_with_flag(handles.F3,'psal','pres',[1:n_prof],'LineWidth',0.5);
			axes(handles.ax2)
		end
	end
end
handles.choice_previous=choice;
guidata(hObject,handles);
%disp('je passe dans Parameter callback3')

% --- Executes during object creation, after setting all properties.
function parameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function profile_Callback(hObject, eventdata, handles)
% hObject    handle to profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of profile as text
%        str2double(get(hObject,'String')) returns contents of profile as a double

vert_sampl=get(handles.vertical_sampling,'String');
choice_vert=get(handles.vertical_sampling,'Value');
%disp('je passe dans Profile callback1')
parameters=get(handles.parameter,'String');
choice=get(handles.parameter,'Value');
%keyboard

if ~isempty(strfind(vert_sampl{choice_vert},'Primary Sampling'))
    vertical_sampling_scheme='Primary sampling';
end
if ~isempty(strfind(vert_sampl{choice_vert},'Near Surface Sampling'))
    vertical_sampling_scheme='Near-surface sampling';
end
k=(get(handles.profile,'String'));
[d,r]=strtok(k,':');
if str2double(d)<1
    k=['[1' r];
end
[d,r]=strtok(k,':');
if isempty(strfind(r,'end'))==0
    r=strrep(r,'end',num2str(length(handles.file_list)));
else
    if str2double(strrep(r,':',''))>length(handles.file_list)
        k=[ d ':end'];
    end
end

%eval(['thelist=handles.file_list(' k ');']);

% if k<1
%     k=1;
% end
% if k>length(handles.file_list)
%     k=length(handles.file_list);
% end

eval(['thelist=handles.file_list(' k ');']);
if strcmp(parameters{choice},'Oxygene')||strcmp(parameters{choice},'Oxygene(Adjusted)')
eval(['thelist=handles.file_listB(' k ');']);
end
eval(['thelistC=handles.file_list(' k ');']);
%keyboard
if length(thelist)>1
    Param='';
    %vertical_sampling_scheme='Primary sampling';
    display('Ploting data...')
    [F,Dim,thelist_ext]=create_multi_from_filelist(handles.floatname,handles.dacname,handles.DIR_FTP,thelist,vertical_sampling_scheme,Param);
	[F3,Dim3,thelist_ext3]=create_multi_from_filelist(handles.floatname,handles.dacname,handles.DIR_FTP,thelistC,vertical_sampling_scheme,Param);
    G='';
    F2='';
    thelist=thelist_ext
else

    k=[str2num(k)-1:str2num(k)+1];
    k(k<1)=1;
    k(k>length(handles.file_list))=length(handles.file_list);
    k=unique(k);
    kstr=[num2str(k(1)) ':' num2str(k(end))];
    eval(['thelist2=handles.file_list(' kstr ');']);
	if strcmp(parameters{choice},'Oxygene')||strcmp(parameters{choice},'Oxygene(Adjusted)')
		eval(['thelist2=handles.file_listB(' kstr ');']);
	end
     Param='';
     
    [F,Dim,G] = read_netcdf_allthefile([handles.DIR_FTP handles.dacname '/' handles.floatname '/profiles/' thelist{1}]);
    [F2,Dim2,thelist_ext2]=create_multi_from_filelist(handles.floatname,handles.dacname,handles.DIR_FTP,thelist2,vertical_sampling_scheme,Param);
	[F3,Dim3,G3] = read_netcdf_allthefile([handles.DIR_FTP handles.dacname '/' handles.floatname '/profiles/' thelistC{1}]);
    is_primary=findstr_tab(F.vertical_sampling_scheme.data,vertical_sampling_scheme);
    [F,Dim] = extract_profile_dim(F,Dim,'N_PROF',is_primary);
	[F3,Dim3] = extract_profile_dim(F3,Dim3,'N_PROF',is_primary);
    F2 = format_flags_char2num(F2);
    F2 = replace_fill_bynan(F2);
	%disp('je passe dans Profile callback2')
	

end

F = format_flags_char2num(F);
F = replace_fill_bynan(F);
F3 = format_flags_char2num(F3);
F3 = replace_fill_bynan(F3);
handles.thelist=thelist;
handles.thelistC=thelistC;
handles.F=F;
handles.F2=F2;
handles.F3=F3;
handles.Dim=Dim;
handles.G=G;
guidata(hObject,handles);
parameter_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function profile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prof_number=get(handles.profile,'String');
if isempty(prof_number)==0
    if isempty(strfind(prof_number,':'))==0||isempty(strfind(prof_number,','))
        [d,r]=strtok(prof_number,':');
        [d,r]=strtok(prof_number,',');
		d=strrep(d,'[','');
        k=str2double(d)-1;
    else
        k=str2num(prof_number)-1;
    end
    if k<1
        k=1;
    end
    if k>length(handles.file_list)
        k=length(handles.file_list);
    end
    prof_number=num2str(k);   set(handles.profile,'String',prof_number);
    profile_Callback(hObject, eventdata, handles)
    
else
    set(handles.profile,'String','1');
    profile_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in NEXT.
function NEXT_Callback(hObject, eventdata, handles)
% hObject    handle to NEXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%keyboard
prof_number=get(handles.profile,'String');
%keyboard
if isempty(prof_number)==0
    
    if isempty(strfind(prof_number,':'))==0||isempty(strfind(prof_number,','))
	    
	    r=prof_number;
		r=strrep(r,']','');
		r=strrep(r,'[','');
		while isempty(strfind(r,':'))==0
        [d,r]=strtok(r,':');
		if isempty(r)==0
		r=r(2:end);
		end
		end
		while isempty(strfind(r,','))==0
        [d,r]=strtok(r,',');
		if isempty(r)==0
		    r=r(2:end);
		end
		end
		
        if isempty(strfind(r,'end'))==0
            r=strrep(r,'end',num2str(length(handles.file_list)));
        end
        k=str2double(r)+1;
    else
        k=str2num(prof_number)+1;
    end
    if k<1
        k=1;
    end
    if k>length(handles.file_list)
        k=length(handles.file_list);
    end
    prof_number=num2str(k);
    set(handles.profile,'String',prof_number);
    profile_Callback(hObject, eventdata, handles)
    
else
    set(handles.profile,'String','1');
    profile_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in validate.
function validate_Callback(hObject, eventdata, handles)
% hObject    handle to validate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save the file
% k=str2double(get(handles.profile,'String'));
% if k<1
%     k=1;
% end
% if k>length(handles.file_list)
%     k=length(handles.file_list);
% end


thedate=datestr(now,'yyyymmddHHMMSS');
parameters=get(handles.parameter,'String');
choice=get(handles.parameter,'Value');
param='';


if strcmp(parameters{choice},'Temperature')
    param='temp';
    thepres='pres';
end
if strcmp(parameters{choice},'Salinity')
    param='psal'; thepres='pres';
end
if strcmp(parameters{choice},'Oxygene')
    param='doxy'; thepres='pres';
end
% if strcmp(parameters{choice},'Temperature(Adjusted)')
%     param='temp_adjusted'; thepres='pres_adjusted';
% end
% if strcmp(parameters{choice},'Salinity(Adjusted)')
%     param='psal_adjusted';thepres='pres_adjusted';
% end
if isempty(param)
    helpdlg(['no validation possible for ' parameters{choice}])
     return
end
vert_sampl=get(handles.vertical_sampling,'String');
choice_vert=get(handles.vertical_sampling,'Value');
%keyboard
% vertical sampling
if ~isempty(strfind(vert_sampl{choice_vert},'Primary Sampling'))
    vertical_sampling_scheme='Primary sampling';
end
if ~isempty(strfind(vert_sampl{choice_vert},'Near Surface Sampling'))
    vertical_sampling_scheme='Near-surface sampling';
end


param_qc=[param '_qc'];
for k=1:length(handles.thelist)
    
    [F,DimF,Glo]= read_netcdf_allthefile([handles.DIR_FTP handles.dacname '/' handles.floatname '/profiles/' handles.thelist{k}]);
	
    [G,DimG]=extract_profile_dim(handles.F,handles.Dim,'N_PROF',k);
    
    % verifie qu'on est bien sur le meme cycle
    if F.cycle_number.data~=G.cycle_number.data|F.direction.data~=G.direction.data
        error('pas le meme cycle')
    end
    
    is_primary=findstr_tab(F.vertical_sampling_scheme.data,vertical_sampling_scheme);
    
    if sum(is_primary)~=1
        error('pas le vertical sampling scheme requis')
    end
    n_prof=find(is_primary==1);
    %keyboard
    sizF=size(F.pres.data(n_prof,:),2);
    G=format_flags_num2char(G);
    G=replace_nan_byfill(G);
    
    old_qc = F.(param_qc).data(n_prof,:);
    new_qc = G.(param_qc).data(1:sizF);
    
    
    %trouve les start_pres,stop_pres et  prev_val
    isneq=find(old_qc~=new_qc);
    if ~isempty(isneq)
        clear start_pres stop_pres prev_val
        start_pres(1)=F.(thepres).data(n_prof,isneq(1));
        stop_pres(1)=F.(thepres).data(n_prof,isneq(end));
        all_qc=str2num(old_qc(isneq)');
        prev_val(1)=all_qc(1);
        
        kk=1;
        for jk=2:length(isneq)
            
            if (isneq(jk)~=isneq(jk-1)+1)| all_qc(jk)~=all_qc(jk-1)
                kk=kk+1;
                stop_pres(kk-1) = F.(thepres).data(n_prof,isneq(jk-1));
                start_pres(kk) = F.(thepres).data(n_prof,isneq(jk));
                stop_pres(kk) =F.(thepres).data(n_prof,isneq(end));
                prev_val(kk) =all_qc(jk);
            end
        end
        
        
        F.(param_qc).data(n_prof,:)=G.(param_qc).data(1:sizF);
        
        % PROFILE_QC
        %keyboard
        [F,GlobQC,ValnumQC]=check_profile_qc(F);
        
        for jk=1:length(start_pres)
            
            % HISTORY fields
            new_hist = DimF.n_history.dimlength+1;
            
            allfields = fieldnames(F);
            
            ii = strfind(allfields,'history_');
            is_history = find(~cellfun('isempty',ii));
            
            F = check_FirstDimArray_is(F,'N_HISTORY');
            
            if DimF.n_history.dimlength~=0
                % remplit un historique de plus avec des FillValue
                [F_ex,DimF_ex]=extract_profile_dim(F,DimF,'N_HISTORY',1);
                %keyboard
                for ik = is_history'
                    oneChamp =allfields{ik};
                    ii=F_ex.(oneChamp).data~=F_ex.(oneChamp).FillValue_;
                    F_ex.(oneChamp).data(ii)=F_ex.(oneChamp).FillValue_;
                end
                
                [F,DimF] = cat_profile_dim(F,F_ex,DimF,DimF_ex,'N_HISTORY');
            else
                for ik = is_history'
                    oneChamp =allfields{ik};
                    siz(1)=1;
                    for tk=2:length(F.(oneChamp).dim)
                        siz(tk) = DimF.(lower(F.(oneChamp).dim{tk})).dimlength;
                    end
                    F.(oneChamp).data = repmat(F_ex.(oneChamp).FillValue_,siz);
                    DimF.n_history.dimlength=1;
                end
            end
            
            
            % remplit l'historique par la valeur appropriee
            step='ARSQ';
            l_in=length(step);
            F.history_step.data(new_hist,n_prof,1:l_in)=step;
            
            institution='IF';
            l_in=length(institution);
            F.history_institution.data(new_hist,n_prof,1:l_in)=institution;
            
            soft='SCOO';
            l_so=length(soft);
            F.history_software.data(new_hist,n_prof,1:l_so)=soft;
            
            soft_release='1.4';
            l_so_r=length(soft_release);
            F.history_software_release.data(new_hist,n_prof,1:l_so_r)=soft_release;
            
            F.history_date.data(new_hist,n_prof,:)=thedate;
            
            action='CF';
            l_ac=length(action);
            F.history_action.data(new_hist,n_prof,1:l_ac)=action;
            
            parameter=upper(param);
            l_pa=length(parameter);
            F.history_parameter.data(new_hist,n_prof,1:l_pa)=parameter;
            
            
            F.history_start_pres.data(new_hist,n_prof)=start_pres(jk);
            F.history_stop_pres.data(new_hist,n_prof)=stop_pres(jk);
            
            
            F.history_previous_value.data(new_hist,n_prof)=prev_val(jk);
        end
        F = check_FirstDimArray_is(F,'N_PROF');
        ficout=[handles.DIR_FTP handles.dacname '/' handles.floatname '/profiles/' handles.thelist{k}];
        
        
        %ficout=['/tmp/' handles.file_list{k}];
        create_netcdf_allthefile(F,DimF,ficout,Glo);
    end
end

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
profile_Callback(hObject, eventdata, handles)


% --- Executes on selection change in vertical_sampling.
function vertical_sampling_Callback(hObject, eventdata, handles)
% hObject    handle to vertical_sampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vertical_sampling contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vertical_sampling
profile_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function vertical_sampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vertical_sampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

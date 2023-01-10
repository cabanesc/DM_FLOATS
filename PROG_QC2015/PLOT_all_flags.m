% -========================================================
%   USAGE : PLOT_all_flags(floatname,dacname,varargin)
%   PURPOSE : plot current flag in netcdf files
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char) e.g.  'coriolis'
%
%   OPTIONNAL INPUT :
%   DATAREP
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created (2016) ccabanes
%
%   CALLED SUBROUTINES: 
% -------------------------------------
% GIT BRANCH: PourOXY
% ========================================================

function PLOT_all_flags(floatname,dacname,varargin)


%close all

C=load_configuration('config.txt');
n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);
ERASE=0;
VPN=1;
if isfield(s,'ERASE')==1;ERASE=s.ERASE;end;
if isfield(s,'VPN')==1;VPN=s.VPN;end;

DIR_FTP_CORIOLIS=C.DIR_FTP_CORIOLIS;
CONFIG.FILE_TOPO=C.FILE_TOPO;
CONFIG.DIR_FTP=C.DIR_FTP;
CONFIG.DIR_PLOT=C.DIR_PLOT;

% FLOTEUR ANALYSE
CONFIG.floatname=floatname;
CONFIG.dacname=dacname;

plotpath=[CONFIG.DIR_PLOT 'preliminaire/' CONFIG.floatname '/'];

if ~exist(plotpath,'dir')
    mkdir(plotpath);
end

[file_list]=select_float_files_on_ftp(CONFIG.floatname,CONFIG.dacname,CONFIG.DIR_FTP,'C');
vertical_sampling_scheme='Primary sampling';
Param.cycle_number.name='CYCLE_NUMBER';
Param.temp.name='TEMP';
Param.temp_qc.name='TEMP_QC';
Param.temp_adjusted.name='TEMP_ADJUSTED';
Param.temp_adjusted_qc.name='TEMP_ADJUSTED_QC';
Param.psal.name='PSAL';
Param.psal_qc.name='PSAL_QC';
Param.psal_adjusted.name='PSAL_ADJUSTED';
Param.psal_adjusted_qc.name='PSAL_ADJUSTED_QC';
Param.psal_adjusted_error.name='PSAL_ADJUSTED_ERROR';
Param.pres.name='PRES';
Param.pres_qc.name='PRES_QC';
Param.pres_adjusted.name='PRES_ADJUSTED';
Param.pres_adjusted_qc.name='PRES_ADJUSTED_QC';
Param.juld.name='JULD';
Param.data_mode.name='DATA_MODE';
[F,Dim,thelist_ext]=create_multi_from_filelist(CONFIG.floatname,CONFIG.dacname,CONFIG.DIR_FTP,file_list,vertical_sampling_scheme,Param);

F=format_flags_char2num(F);
fflag=plot_one_flag(F,'psal')
eval(['print(fflag, ''-dpng'',''' plotpath CONFIG.floatname  '_current_flags_PSAL.png'')'])

fflag=plot_one_flag(F,'temp')
eval(['print(fflag, ''-dpng'',''' plotpath CONFIG.floatname  '_current_flags_TEMP.png'')'])

fflag=plot_one_flag(F,'pres')
eval(['print(fflag, ''-dpng'',''' plotpath CONFIG.floatname  '_current_flags_PRES.png'')'])

  
F = replace_fill_bynan(F);
diff=F.pres_adjusted.data-F.pres.data;
SALINI=F.psal.data;

if meanoutnan(meanoutnan(diff))~=0&~isnan(meanoutnan(meanoutnan(diff)))

    ij=find(diff~=0&~isnan(diff));
    % compute conductivity from the salinity and raw pressure
    cndr = sw_cndr(F.psal.data(ij),F.temp.data(ij),F.pres.data(ij));
    PRESINI=F.pres_adjusted.data;
    
    
    % recompute the salinity from conductivity and adjusted pressure
    sal = sw_salt(cndr,F.temp.data(ij),F.pres_adjusted.data(ij));
    sal(isnan(sal))=NaN;
    SALINI(ij)=sal;
end
F.psal.data=SALINI;
N_PROF_this_file = size(F.cycle_number.data,1);
param='psal';
paramad='psal_adjusted';
paramqc='psal_qc';
paramadqc='psal_adjusted_qc';

mean_diff = NaN*zeros(N_PROF_this_file,1);

%isnanprof=isnan(F.(param).data)|isnan(F.(paramad).data)|F.(param).data>9999999|F.(paramad).data>9999999;
isnanprof=isnan(F.(param).data)|isnan(F.(paramad).data);
F.(paramad).data(isnanprof)=NaN;
F.(param).data(isnanprof)=NaN;

%F.(param).data(F.(param).data>9999999)=NaN;
%F.(paramad).data(F.(paramad).data>9999999)=NaN;

mean_diff = meanoutnan(F.(paramad).data,2)-meanoutnan(F.(param).data,2);
mean_error= meanoutnan(F.([paramad '_error']).data,2);

b=double([F.psal_adjusted.data-F.psal.data;F.psal_adjusted.data(end,:)-F.psal.data(end,:)]);
if exist(plotpath)==0
mkdir(plotpath)
end

figure
subplot(2,1,1)
hold on
box on
grid on
set(gca,'Fontsize',12)
set(gca,'FontWeight','bold')
%plot(mean_diff)
%keyboard
%pcolor(a,[1:size(F.psal_adjusted_qc.data,2)],b')
%set(gca,'YDir','reverse')
%shading('flat')
%ylabel('Mean PSAL Correction','Fontsize',10)
xlabel('Cycle Number','interpreter','none','Fontsize',10)

%keyboard
plot(F.cycle_number.data,mean_diff,'m+','LineWidth',4)
yax=max(abs(get(gca,'YLim')));
axis([0 F.cycle_number.data(end)+1+0.5 -yax yax])
ylabel('psu','interpreter', 'none','Fontsize',10)
xlabel('Cycle Number','interpreter','none','Fontsize',10)
title([{'Vertical Mean PSAL correction'},{'(psal_adjusted - psal)'}],'interpreter', 'none')

% subplot(2,1,2)
% hold on
% box on
% grid on
% set(gca,'Fontsize',12)
% set(gca,'FontWeight','bold')
% %plot(mean_diff)
% %keyboard
% %pcolor(a,[1:size(F.psal_adjusted_qc.data,2)],b')
% %set(gca,'YDir','reverse')
% %shading('flat')
% %ylabel('Mean PSAL Correction','Fontsize',10)
% xlabel('Cycle Number','interpreter','none','Fontsize',10)
% 
% %keyboard
% plot(F.cycle_number.data,mean_error,'b-+','LineWidth',4)
% 
% ylim=get(gca,'YLim');
% %yax=max(abs(get(gca,'YLim')));
% axis([0 F.cycle_number.data(end)+1+0.5 ylim(1) ylim(2)])
% ylabel('psu','interpreter', 'none','Fontsize',10)
% xlabel('Cycle Number','interpreter','none','Fontsize',10)
% title(['Mean PSAL_ADJUSTED ERROR'],'interpreter', 'none')
%c=colorbar('Location', 'SouthOutside')
%xlabel(c,'DM PSAL correction (PSU)')
plotpath=[CONFIG.DIR_PLOT 'preliminaire/' CONFIG.floatname '/'];
if exist(plotpath)==0
mkdir(plotpath)
end
%caxis([-0.03 0])

eval(['print -dpng ' plotpath CONFIG.floatname  '_current_adjustement.png'])



% compute the drift rate over each 10 cycles
%------------------------
disp('')
disp('DRIFT RATE over 10 cycles')
disp('-------------------------')
tt=[1:10:length(mean_diff)];
c=F.cycle_number.data;
b=F.juld.data;
if length(tt)>1
	for i=2:length(tt)

		 r=-(mean_diff(tt(i))-mean_diff(tt(i-1)))./(b(tt(i))-b(tt(i-1)))*365;
		 disp(['cycles ' num2str(c(tt(i-1))) ' to ' num2str(c(tt(i))) ': ' num2str(r) ' PSU/yr'])
	end
end


% last cycle
%-----------
disp('-------------------------')

disp('')
disp(['LAST CYCLE: ' num2str(F.cycle_number.data(end)) ' - ' datestr(datenum('19500101','yyyymmdd')+F.juld.data(end))])

v=find(findstr_tab(cellstr(F.data_mode.data),'D'));
if isempty(v)==0
disp(['LAST CYCLE in DM: ' num2str(F.cycle_number.data(v(end))) ])
end

% thedate=input('enter date yyyymmdd: ','s');
% thedatenum=datenum(thedate,'yyyymmdd');
% thedatefloat=datenum('19500101','yyyymmdd')+F.juld.data;
% ff=find(thedatefloat>=thedatenum);
% if isempty(ff)==0
% disp(['Cycle number: ' num2str(F.cycle_number.data(ff(1))) ])
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figure1=plot_one_flag(F,param)

paramqc=[param '_qc'];
paramad=[param '_adjusted'];
paramadqc=[param '_adjusted_qc'];

figure1=figure;
subplot(2,1,1)
set(gca,'Fontsize',12)
set(gca,'FontWeight','bold')
hold on
box on
w(:,:)=[1,1,1;0,0,1;0.3,1,0.3;1,1,0.3;1,0.7,0.3;1,0.3,0.3];
F.(paramqc).data(F.(paramqc).data==999)=-1;
a=[F.cycle_number.data;F.cycle_number.data(end)+1];
b=double([F.(paramqc).data;F.(paramqc).data(end,:)]);
bn=b;bn(bn==-1)=NaN;
d3= find( meanoutnan(bn,2)==3);
d4= find( meanoutnan(bn,2)==4);



if isempty(d3)==0
disp(upper(param))    
disp(['flag3: from cycle ' num2str(a(d3(1)))])
end
if isempty(d4)==0
disp(upper(param))      
disp(['flag4: from cycle ' num2str(a(d4(1)))])
end

A=repmat(a',[1,size(b,2)]);
c=[1:size(F.(paramqc).data,2)];
C=repmat(c, [size(b,1),1]);
ii=find(b==0);
scatter(A(ii),C(ii),10,'*b')
ii=find(b==1);
scatter(A(ii),C(ii),10,'*g')
ii=find(b==2);
scatter(A(ii),C(ii),10,'*y')
ii=find(b==3);
scatter(A(ii),C(ii),10,'*m')
ii=find(b==4);
scatter(A(ii),C(ii),10,'*r')
set(gca,'YDir','reverse')
ylabel('level')
%xlabel('Cycle Number','interpreter','none')
title([upper(paramqc) ' in the netcdf file'],'interpreter', 'none')

axis([0 F.cycle_number.data(end)+1+0.5 0 size(F.(paramqc).data,2)])
grid on


subplot(2,1,2)
set(gca,'Fontsize',12)
set(gca,'FontWeight','bold')
hold on
box on
w(:,:)=[1,1,1;0,0,1;0.3,1,0.3;1,1,0.3;1,0.7,0.3;1,0.3,0.3];
F.(paramqc).data(F.(paramqc).data==999)=-1;
a=[F.cycle_number.data;F.cycle_number.data(end)+1];
b=double([F.(paramadqc).data;F.(paramadqc).data(end,:)]);
bn=b;bn(bn==999)=NaN;
d3= find( meanoutnan(bn,2)==3);
d4= find( meanoutnan(bn,2)==4);


if isempty(d3)==0
disp(upper(paramad))   
disp(['flag3: from cycle ' num2str(a(d3(1)))])
end
if isempty(d4)==0
disp(upper(paramad))      
disp(['flag4: from cycle ' num2str(a(d4(1)))])
end

A=repmat(a',[1,size(b,2)]);
c=[1:size(F.(paramadqc).data,2)];
C=repmat(c, [size(b,1),1]);
%pcolor([1:size(CHECK.(param).adjqc,1)],[1:size(CHECK.(param).adjqc,2)],double(CHECK.(param).adjqc'))
%pcolor(a,[1:size(CHECK.(param).qc,2)],b')
ii=find(b==0);
scatter(A(ii),C(ii),10,'*b')
ii=find(b==1);
scatter(A(ii),C(ii),10,'*g')
ii=find(b==2);
scatter(A(ii),C(ii),10,'*y')
ii=find(b==3);
scatter(A(ii),C(ii),10,'*m')
ii=find(b==4);
scatter(A(ii),C(ii),10,'*r')
set(gca,'YDir','reverse')
%shading('flat')
ylabel('level')
xlabel('Cycle Number','interpreter','none')
title([upper(paramadqc) ' in the netcdf file'],'interpreter', 'none')
axis([0 F.cycle_number.data(end)+1+0.5 0 size(F.(paramadqc).data,2)])

grid on

X=get(gca,'XLIM');
% Create rectangle
annotation(figure1,'textbox',...
    [0.926 0.60 0.056 0.05],...
    'string','Flags','FontWeight','bold','Color',[0 0 0],'HorizontalAlignment','center','EdgeColor','none');

annotation(figure1,'textbox',...
    [0.926 0.55 0.056 0.05],...
    'BackgroundColor',[1 0 0],'string','4','FontWeight','bold','Color',[0.6 0.6 0.6],'HorizontalAlignment','center');

annotation(figure1,'textbox',...
    [0.926 0.50 0.056 0.05],...
    'BackgroundColor',[1 0 1],'string','3','FontWeight','bold','Color',[0.6 0.6 0.6],'HorizontalAlignment','center');

annotation(figure1,'textbox',...
    [0.926 0.45 0.056 0.05],...
    'BackgroundColor',[1 1 0],'string','2','FontWeight','bold','Color',[0.6 0.6 0.6],'HorizontalAlignment','center');

annotation(figure1,'textbox',...
    [0.926 0.40 0.056 0.05],...
    'BackgroundColor',[0 1 0],'string','1','FontWeight','bold','Color',[0.6 0.6 0.6],'HorizontalAlignment','center');

annotation(figure1,'textbox',...
    [0.926 0.35 0.056 0.05],...
    'BackgroundColor',[0 0 1],'string','0','FontWeight','bold','Color',[0.6 0.6 0.6],'HorizontalAlignment','center');




function plot_float_tpot(floatname,dacname,CONFIG,varargin)
% -========================================================
%   USAGE : plot_float_diag(floatname,dacname,CONFIG,varargin)
%   PURPOSE : trace les données d'un flotteur sur une carte a un niveau theta donne
% -----------------------------------
%   INPUT :
%     IN1   (class)  -comments-
%             additional description
%     IN2   (class)  -comments-
%
%   OPTIONNAL INPUT :
%    OPTION1  (class)  -comments-
% -----------------------------------
%   OUTPUT :
%     OUT1   (class)  -comments-
%             additional description
%     OUT2   (class)  -comments-
%             additional description
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx 
%   CALLED SUBROUTINES: none
% ========================================================
n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);
accepted_options={'MarkerSize','MarkerFaceColor','MarkerSpec'};
accepted_Param={'','psal_mean','thedates','pres_mean','longitude','latitude'};

thefield = fieldnames(s);
for k=1:length(thefield)
    if ismember(thefield{k},accepted_options)==0
    error(['Unreconized option ' thefield{k}])
    end
end
if isfield(s,'MarkerSize')==0
s.MarkerSize=20; %default
end
if isfield(s,'MarkerSpec')==0
s.MarkerSpec='o'; %default
end
if isfield(s,'Isbest')==0
    s.Isbest=0; %default
end



hold on
box on
grid on


ModeClim=get(gca,'CLimMode');
ht=get(gca,'Title');
thetitle=get(ht,'string');
if s.Isbest==0
thestr='(raw)';
elseif s.Isbest==1
thestr='(best)';
elseif s.Isbest==2
thestr='(adj)';
end
if ~isempty(thetitle)&isempty(strfind(thetitle,'float data'))
title([thetitle ' and ' floatname ' float data' thestr])
end
vertical_sampling_scheme='Primary sampling';
IncludeDescProf=1;
[file_list] = select_float_files_on_ftp(floatname,dacname,CONFIG.DIR_FTP,'C',IncludeDescProf);
[F,Dim,thelist_ext2]=create_multi_from_filelist(floatname,dacname,CONFIG.DIR_FTP,file_list,vertical_sampling_scheme,'');

% filename=[CONFIG.DIR_FTP  dacname '/' floatname '/' floatname '_prof.nc'];
% F=read_netcdf_allthefile(filename);
F = replace_fill_bynan(F);
F = format_flags_char2num(F);
F.psal.data(F.psal_qc.data>3)=NaN;
if s.Isbest==1
	F =construct_best_param(F ,{'temp','pres','psal'},F);
	F.psal=F.psal_best;F.temp=F.temp_best;F.pres=F.pres_best;
	F.psal.data(F.psal_qc.data>1)=NaN;
end
if s.Isbest==2
	%F =construct_best_param(F ,{'temp','pres','psal'},F);
	F.psal=F.psal_adjusted;F.temp=F.temp_adjusted;F.pres=F.pres_adjusted;
	F.psal.data(F.psal_adjusted_qc.data>3)=NaN;
end
F.tpot.data = sw_ptmp(F.psal.data,F.temp.data,F.pres.data,0);
POUB=[];
 if CONFIG.ONTHETA==1
[POUB,F] = find_psal_on_theta(F, CONFIG.TPOT_MIN, CONFIG.TPOT_MAX, CONFIG.DEPTH_MIN,POUB);
else
[POUB,F] = find_psal_on_theta(F, CONFIG.P_MIN, CONFIG.P_MAX,CONFIG.DEPTH_MIN,POUB);
end
thedate = datevec((F.juld.data+datenum('19500101','yyyymmdd')));
siz=size(thedate,1);
ll=[thedate(:,1),ones(siz,2),zeros(siz,3)];
F.thedates.data = thedate(:,1)+etime(thedate,ll)./(3600*24*365.25);
F.tpot_min.data=min(F.tpot.data');
F.pres_max.data=max(F.pres.data');

ip=find(F.latitude.data>=CONFIG.VEC_REG(3)&F.latitude.data<=CONFIG.VEC_REG(4)&F.longitude.data<=CONFIG.VEC_REG(2)&F.longitude.data>=CONFIG.VEC_REG(1));
for il=1:length(ip)
plot(F.psal.data(ip(il),:)',F.tpot.data(ip(il),:)',s.MarkerFaceColor);
end   


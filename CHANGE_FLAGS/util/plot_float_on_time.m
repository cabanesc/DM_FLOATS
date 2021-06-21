function plot_float_on_time(floatname,dacname,CONFIG,Param,varargin)
% -========================================================
%   USAGE : plot_float_on_time(floatname,dacname,CONFIG,Param,varargin)
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
if ismember(Param,accepted_Param)==0
   error(['Unreconized Param ' Param])
end
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
if isempty(Param)==1&~isfield(s,'MarkerFaceColor')
   s.MarkerFaceColor='m'; % default
end
hold on
box on
grid on


ModeClim=get(gca,'CLimMode')
ht=get(gca,'Title');
thetitle=get(ht,'string')
if ~isempty(thetitle)&isempty(strfind(thetitle,'float data'))
title([thetitle ' and ' floatname ' float data'])
end

filename=[CONFIG.DIR_FTP  dacname '/' floatname '/' floatname '_prof.nc'];
F=read_netcdf_allthefile(filename);
F = replace_fill_bynan(F);
F = format_flags_char2num(F);
F.psal.data(F.psal_qc.data>2)=NaN;
F.tpot.data = sw_ptmp(F.psal.data,F.temp.data,F.pres.data,0);
POUB=[];
 if CONFIG.OnTheta==1
[POUB,F] = find_psal_on_theta(F, CONFIG.TPOT_MIN, CONFIG.TPOT_MAX, CONFIG.MIN_DEPTH,POUB);
else
[POUB,F] = find_psal_on_theta(F, CONFIG.P_MIN, CONFIG.P_MAX,POUB);
end
thedate = datevec((F.juld.data+datenum('19500101','yyyymmdd')));
siz=size(thedate,1);
ll=[thedate(:,1),ones(siz,2),zeros(siz,3)];
F.thedates.data = thedate(:,1)+etime(thedate,ll)./(3600*24*365.25);
F.tpot_min.data=min(F.tpot.data');
F.pres_max.data=max(F.pres.data');
   
plot(F.thedates.data,F.psal_mean.data,'-m','LineWidth',1)
scatter(F.thedates.data,F.psal_mean.data,s.MarkerSize,F.(Param).data,s.MarkerSpec,'filled')
scatter(F.thedates.data,F.psal_mean.data,s.MarkerSize,['m' s.MarkerSpec])


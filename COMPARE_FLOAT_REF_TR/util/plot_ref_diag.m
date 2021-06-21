function plot_ref_diag(DATA_BASE,CONFIG,varargin)
% -========================================================
%   USAGE : plot_ref_diag(DATA_BASE,CONFIG,varargin)
%   PURPOSE : trace les diagramme theta/S des donnees de reference
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
%  if ismember(Param,accepted_Param)==0
%     error(['Unreconized Param ' Param])
%  end
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
%  if ~isfield(s,'MarkerFaceColor')
%     s.MarkerFaceColor='m'; 
%  end

if CONFIG.CHECK_REF==1
switch DATA_BASE
    case 'CTD'
	filelog=[CONFIG.DIR_PLOT '/log/log_bad_ref_ctd_' CONFIG.VERSION_OWC '.txt'];

    case 'ARGO' 
	filelog=[CONFIG.DIR_PLOT '/log/log_bad_ref_argo_' CONFIG.VERSION_OWC '.txt'];

end
fw=fopen(filelog,'a');
end
ModeClim=get(gca,'CLimMode');
ht=get(gca,'Title');
thetitle=get(ht,'string');



%keyboard
if ~isempty(thetitle)
 if isempty(strfind(thetitle,DATA_BASE))
        title(['Theta/S diagram from ARGO and CTD reference databases']);
 else
        %cla(gca)
        ModeClim='auto';
        title(['Theta/S diagram from ' DATA_BASE ' reference database']);
 end
else
title(['Theta/S diagram from ' DATA_BASE ' reference database']);
end

hold on
box on
grid on

region.lonmin=CONFIG.VEC_REG(1);
region.lonmax=CONFIG.VEC_REG(2);
region.latmin=CONFIG.VEC_REG(3);
region.latmax=CONFIG.VEC_REG(4);
DATA_BASE
% Identification des boites wmo à vérifier
switch DATA_BASE
    case 'CTD'
        load([CONFIG.pathwmobox 'constants/wmo_boxes_ctd.mat']);
        wmo_square = find_wmoboxes( region, la_wmo_boxes);
        la_wmo_num=la_wmo_boxes(ismember(la_wmo_boxes(:,1),wmo_square)&la_wmo_boxes(:,2)==1,:);
        rep='climatology/historical_ctd/ctd_';
    case 'ARGO'
        load([CONFIG.pathwmobox 'constants/wmo_boxes_argo.mat']);
        wmo_square = find_wmoboxes( region, la_wmo_boxes);
        la_wmo_num=la_wmo_boxes(ismember(la_wmo_boxes(:,1),wmo_square)&la_wmo_boxes(:,4)==1,:);
        rep='climatology/argo_profiles/argo_';
end

[nbox,null]=size(la_wmo_num);

if strcmp(CONFIG.VEC_PSAL_AUTO,'on')==0
caxis(CONFIG.VEC_PSAL)
end


if isfield(CONFIG,'YEAR_MIN')==0
CONFIG.YEAR_MIN=-9999;
end
if isfield(CONFIG,'YEAR_MAX')==0
CONFIG.YEAR_MAX=9999;
end
alldates=[];
switch ModeClim
case 'auto'
themin=999999;
themax=-999999;
case 'manual'
color_lim=get(gca,'Clim');
themin=color_lim(1);
themax=color_lim(2);
end

hold on
box on
grid on
cmpa=summer(CONFIG.YEAR_MAX-CONFIG.YEAR_MIN+1);
colormap(cmpa)
theyear=[CONFIG.YEAR_MIN:CONFIG.YEAR_MAX];
for ibox=1:nbox % boucle sur toutes les boites
    clear CTD
    PLOT=[];
    wmo_file= [CONFIG.pathwmobox rep num2str(la_wmo_num(ibox,1))];
    [CTD] = read_wmo(wmo_file,CONFIG.ctd_to_exclude,DATA_BASE);
     CTD.tpot_min.data=min(CTD.tpot.data');
     CTD.pres_max.data=max(CTD.pres.data');
    [PLOT,CTD] = find_psal_on_theta(CTD, CONFIG.TPOT_MIN, CONFIG.TPOT_MAX, CONFIG.DEPTH_MIN,PLOT);
	
	ideC1=~isnan(CTD.psal.data);
	stdbase1=std(CTD.psal.data(ideC1));
	meanbase1=mean(CTD.psal.data(ideC1));
	isout1=abs(CTD.psal.data-meanbase1)>10*stdbase1;
		
		
	if CONFIG.CHECK_REF==1
		
		% log des problemes dans un fichiers log
		findout=find(sum(isout1,2)>=1);
		if isempty (findout)==0
		disp('')
		disp('WARNING')
		disp([num2str(length(findout)) ' Salinity outlier(s) are found in the reference database'])
		disp(['see: ' filelog])
		disp('')
		
		end
		
		for ilog=1:length(findout)
		[puob,ipsal_bad_value]=max(abs(CTD.psal.data(findout(ilog),:)-meanbase1));
		
		fprintf(fw,'%s\n',[num2str(CTD.boxes.data(findout(ilog))) ', ' num2str(CTD.profil_orig.data(findout(ilog))) ', ' num2str(CTD.psal.data(findout(ilog),ipsal_bad_value)) ]);
	    end
		% ouliers are NaN;
		CTD.psal.data(isout1)=NaN;
	end

	
    CTD.thedates.data=((CTD.dates.data-datenum('19700101','yyyymmdd'))/365.25+1970);
    ip=find(CTD.latitude.data>=CONFIG.VEC_REG(3)&CTD.latitude.data<=CONFIG.VEC_REG(4)&CTD.longitude.data<=CONFIG.VEC_REG(2)&CTD.longitude.data>=CONFIG.VEC_REG(1)&CTD.thedates.data>=CONFIG.YEAR_MIN&CTD.thedates.data<CONFIG.YEAR_MAX);
    
    for k=1:length(ip)
        
        idate=find(theyear==round(CTD.thedates.data(ip(k),:)));
        thecolomark=cmpa(idate,:);
        %scatter(CTD.psal.data(ip(k),:)',CTD.tpot.data(ip(k),:)',10,repmat(CTD.longitude.data(ip(k)),[1,size(CTD.psal.data(ip(k),:),2)]),'o','filled')
        if isfield(s,'MarkerFaceColor')
         plot(CTD.psal.data(ip(k),:)',CTD.tpot.data(ip(k),:)','Color',s.MarkerFaceColor);
        else
        plot(CTD.psal.data(ip(k),:)',CTD.tpot.data(ip(k),:)','+','Color',thecolomark);
        end
      
    end
    if isempty(ip)==0
      themin=min(themin,min(CTD.psal.data(ip)));
      themax=max(themax,max(CTD.psal.data(ip)));
    end
end 
themin
themax
%keyboard
  caxis([themin,themax])
%  colorbar
%  a=get(gca,'XLim');
%  b=get(gca,'YLim');
%  text(a(2)+(a(2)-a(1))/4.7,b(1)+2/3*(b(2)-b(1)),upper(Param),'Rotation',-90,'FontSize',11)
if CONFIG.CHECK_REF==1
fclose(fw)
end
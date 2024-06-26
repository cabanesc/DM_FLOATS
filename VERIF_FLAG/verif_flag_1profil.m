% -========================================================
%   USAGE :   verif_flag_1profil(floatname,dacname,numcycle,varargin)
%   PURPOSE : plot data from a given float
% -----------------------------------
%   INPUT :
%    floatname  (char or cell of char)  e.g. '690258' or {'690258', '3901954'}
%    dacname    (char or cell of char)  e.g.  'coriolis' or {'coriolis', 'bodc'}
%    numcycle   (float array or cell of float array)  e.g. 1 (first cycle 001)  or {1,[1,2]};
%                              see option 'DIRECTION' to consider descending profiles
%    config_campaign_file (char)        e.g. 'config_rrex17.txt'
%   OPTIONNAL INPUT :
%    'CONFIG_OW'  (float)  give the config number for OW and retrieve
%    the value of MAP_USE_PV set in the config file. This value will be used to select the reference data
%    (PARAM.USE_PV=MAP_USE_PV). By default PARAM.USE_PV=1;
%    
%    'DIRECTION' (char)  'A'(default) or 'D' if you want to use ascending or descending profiles
%    'REFERENCE' char)   'ctd'  if reference CTD are used for comparison (default)
%                        'argo' if reference argo data are used for comparison
%    'REXTEND'   (float)  e.g 50: Number of neighboring profiles selected (default: 50)
%    'DATATYPE'  (char)  'raw' if raw PARAM are plotted; 'adj' if PARAM_ADJ are plotted  (default is 'raw')
%    'DATAREP'   (char)   default 'DIR_FTP'  directory where nc files are looked for
%    'DEPHT_ZOOM' (float) 1000 (default) zoom below DEPHT_ZOOM for plotting sig0 and theta/S plots
%    'CHECK_REF'  (logical) 1         : a basic qc is done on the reference data: outliers std*10 are NaN, and a log file with alert is created
%                           0(default): no basic qc (all data are plotted)
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  :  created: V. Thierry - N. David - Juillet 2007
%               revised: C. Lagadec
%               revised: C. Cabanes - 2016
%   CALLED SUBROUTINES:
% ========================================================
function verif_flag_1profil(floatname,dacname, numcycle, varargin)
close all

%Commented by T. Reynaud 08.09.2020
%init_path

if iscell(floatname)==0
floatname=cellstr(floatname);
end
if iscell(dacname)==0
dacname=cellstr(dacname)
end
if iscell(numcycle)==0
numcycle={numcycle};
end

if length(floatname)~=length(dacname)|length(floatname)~=length(numcycle)
    error('floatname/dacname and or numcycle do not have the same length')
end


CONFIG=load_configuration('config.txt');



n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);


% default CONFIG
PARAM.DIRECTION='A';
PARAM.REFERENCE='ctd';
PARAM.REXTEND=50;
PARAM.DATATYPE='raw';
PARAM.DEPTH_ZOOM=1000;
PARAM.CHECK_REF=0;
PARAM.DATAREP='DIR_FTP';


% Input CONFIG
if isfield(s,'DIRECTION')==1;PARAM.DIRECTION=s.DIRECTION;end;
if isfield(s,'REFERENCE')==1;PARAM.REFERENCE=s.REFERENCE;end;
if isfield(s,'DATATYPE')==1;PARAM.DATATYPE=s.DATATYPE;end;
if isfield(s,'REXTEND')==1;PARAM.REXTEND=s.REXTEND;end;
if isfield(s,'DEPTH_ZOOM')==1;PARAM.DEPTH_ZOOM=s.DEPTH_ZOOM;end;
if isfield(s,'CHECK_REF')==1;PARAM.CHECK_REF=s.CHECK_REF;end;
if isfield(s,'DATAREP')==1;PARAM.DATAREP=s.DATAREP;end;
if strcmp(PARAM.DATAREP,'DIR_DM_FILES')
thepostfix='_adj';
else
thepostfix='';
end

if isfield(s,'CONFIG_OW')==1
    if exist([CONFIG.DIR_CODES 'LPO_CODES_ATLN_NEW/ow_config/ow_config_' num2str(s.CONFIG_OW) '.txt'],'file')
        lo_system_configuration = load_configuration( [CONFIG.DIR_CODES 'LPO_CODES_ATLN_NEW/ow_config/ow_config_' num2str(s.CONFIG_OW) '.txt']);
        PARAM.USE_PV= str2double(lo_system_configuration.MAP_USE_PV);
    else
        PARAM.USE_PV=1;
    end
else
    PARAM.USE_PV=1;
end

for ifloat=1:length(floatname)
    
	floatname_i=strtrim(floatname{ifloat});
	dacname_i =strtrim(dacname{ifloat});
	tab=numcycle{ifloat};

	plotpathini=[CONFIG.DIR_PLOT '/verif_flag/'];

	if exist([plotpathini floatname_i]) == 0
		%eval(['!mkdir ' plotpathini floatname_i]);
		mkdir( [plotpathini floatname_i])
	end

		 
	[ntab,ncol]=size(tab);

	for iprf=1:ncol
		wmonum=str2num(floatname_i);
		profnum=tab(iprf);
		plotpath=[plotpathini '/' floatname_i '/'];
		disp(['Flotteur ' floatname_i  '- cycle ' num2str(tab(iprf))]);
		if strcmp(PARAM.REFERENCE,'ctd') == 1
			disp(['Comparison with reference CTD close in space (' num2str(PARAM.REXTEND) ' nearest profiles)']);
			[msg,hf,h1,h2,h3]=cmp_prf_argo_ref(wmonum,dacname_i,profnum,PARAM,CONFIG);
			if strcmp(msg,'ok') == 1
				
				figure(hf)
				titplt=[floatname_i '_prof' num2str(profnum) '_1' thepostfix '.png'];
				eval(['print -dpng ' plotpath titplt]);

				figure(h1)
				titplt=[floatname_i '_prof' num2str(profnum) '_2' thepostfix '.png'];
				eval(['print -dpng ' plotpath 'verif_flag_' titplt]);

				figure(h2)
				titplt=[floatname_i '_prof' num2str(profnum) '_3' thepostfix '.png'];
				eval(['print -dpng ' plotpath 'verif_flag_' titplt]);

				figure(h3)
				titplt=[floatname_i '_prof' num2str(profnum) '_4' thepostfix '.png'];
				eval(['print -dpng ' plotpath 'verif_flag_' titplt]);


			end
		end

		if strcmp(PARAM.REFERENCE,'argo') == 1
			disp(['Comparison with reference Argo close in space (' num2str(PARAM.REXTEND) ' nearest profiles)']);

			[msg,hf,h1,h2,h3]=cmp_prf_argo_ref(wmonum,dacname_i,profnum,PARAM,CONFIG);
			if strcmp(msg,'ok') == 1
			  
				figure(hf)
                set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
				titplt=[floatname_i '_cmpARGO_prof' num2str(profnum) '_1' thepostfix '.png'];
				eval(['print -dpng ' plotpath titplt]);

				figure(h1)
                set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
				titplt=[floatname_i '_cmpARGO_prof' num2str(profnum) '_2' thepostfix '.png'];
				eval(['print -dpng ' plotpath 'verif_flag_' titplt]);

				figure(h2)
                set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
				titplt=[floatname_i '_cmpARGO_prof' num2str(profnum) '_3' thepostfix '.png'];
				eval(['print -dpng ' plotpath 'verif_flag_' titplt]);

				figure(h3)
                set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
				titplt=[floatname_i '_cmpARGO_prof' num2str(profnum) '_4' thepostfix '.png'];
				eval(['print -dpng ' plotpath 'verif_flag_' titplt]);
				
			end
		end
	end

	% conversion of .eps to .pdf
	plotfile=dir([plotpathini '/' floatname_i '/*.eps']);
	fid10=fopen([plotpath '/convert.sh'],'w');
	for k=1:length(plotfile)
		thefile=[plotpathini '/' floatname_i '/' plotfile(k).name];
		thestring=['ps2pdf ' thefile ' ' strrep(thefile,'.eps','.pdf')];
		fprintf (fid10,'%s \n', thestring);
	end
	fclose(fid10);
	eval(['!chmod 777  ' plotpath '/convert.sh'])   
end

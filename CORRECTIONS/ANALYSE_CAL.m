% ========================================================
%   USAGE :  ANALYSE_CAL(flt_name,dacname,varargin)
%
%   PURPOSE : Analyse du fichier CAL (determination de l'offset et de la pente sur differents segments)
%
% -----------------------------------
%   INPUT :
%     flt_name  (char)  -- wmo float name  e.g.: '4900139'
%     dacname   (char)  --   float dac name
%%   OPTIONNAL INPUT :
%    NUMCONFIG (char) : configuration number for OWC (e.g. '149') Needed only if the configuration file is the generic one. (config.txt)
% -----------------------------------
%   OUTPUT :
%   netcdf core D files (single cycle files) that can be distributed on gdac.
% -----------------------------------
%   HISTORY
%  ??/???? : C.Coatanoan: - creation (fichier_nc_dmqc_ovide.m)
%  from 12/2014 : C.Cabanes : -interactive input, use of +libargo, format 3.1
%
%  tested with matlab  8.3.0.532 (R2014a)
%
%  EXTERNAL LIB
%  package +libargo : addpath('dm_float/lib/')
%  seawater  : addpath('dm_float/lib/seawater_330_its90/')
%
%  CONFIGURATION file: config.txt;
%==================================================
function ANALYSE_CAL(flt_name,dacname,varargin)
% INPUT PARAMETERS
n=length(varargin);
if n/2~=floor(n/2)
    error('check the imput arguments')
end
f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);


disp('%%%%%%%')
disp('Delayed-Mode NetCDF files')
% if run in create_dm_files directory
% Commented by T. Reynaud 02.10.2020
% addpath('lib/')
% addpath('lib/seawater_330_its90/')

% INPUT/OUTPUT files

if ischar(flt_name)==0
    flotteur = num2str(flt_name);
else
    flotteur = strtrim(flt_name);
end


if exist(['./paramlog/config_' flotteur '.txt'])
    disp(['CONFIGURATION FILE USED: ./paramlog/config_' flotteur '.txt'])
    
    C=load_configuration(['./paramlog/config_' flotteur '.txt']);
    DIR_FTP= C.DIR_FTP;    % input files directory
    DIR_OW = C.DIR_DATA;     % calibration files directory(cal_$flt_name$.mat files)
else
    disp(['CONFIGURATION FILE USED: config.txt'])
    % default
    numconfig='149';
    if isfield(s,'NUMCONFIG')==1;numconfig=s.NUMCONFIG;end;
    disp(['numconfig =' numconfig])
    C=load_configuration(['../config.txt']);
    %eval(['!cp ../config.txt ./paramlog/config_' flotteur '.txt'])
    DIR_FTP= [C.DIR_FTP dacname '/'];    % input files directory
    DIR_OW = [C.DIR_DATA 'float_calib/CONFIG' numconfig '/' ];     % calibration files directory(cal_$flt_name$.mat files)

end



DIR_CAL=[DIR_OW '/'  ];

C_FILE = load([DIR_CAL '/cal_' num2str(flotteur) '.mat']);
S_FILE = load([DIR_CAL '../../float_source/' num2str(flotteur) '.mat']);

% copie des fichiers cal sur le DIR_CAL


n=size(S_FILE.SAL,2);
Soffset=C_FILE.cal_SAL-S_FILE.SAL;
avg_Soffset=NaN.*ones(1,n);
avg_Soffset_err=NaN.*ones(1,n);
for i=1:n
                ii=[];
                ii=find(isnan(Soffset(:,i))==0);
                if ~isempty(ii)
                    avg_Soffset(i)=mean(Soffset(ii,i));
                    avg_Soffset_err(i)=mean(C_FILE.cal_SAL_err(ii,i));
                else
                    avg_Soffset(i) = NaN;
                    avg_Soffset_err(i) = NaN;
                end
end
conduc= avg_Soffset(~isnan(avg_Soffset));
conduc_err= avg_Soffset_err(~isnan(avg_Soffset_err));

prof_no= C_FILE.PROFILE_NO(~isnan(C_FILE.pcond_factor));
derive_cond= (round((conduc(2:end)-conduc(1:end-1))*1000)/1000)./(prof_no(2:end)-prof_no(1:end-1));


i=1;
segment=[1];
u=2;
while length(segment)<length(derive_cond)-1
    if derive_cond(u)==derive_cond(u-1)%| derive_cond(u)==derive_cond(u+1)
        segment=[segment i];
    else
        i=i+1;
        segment=[segment i];
    end
    u=u+1;
end
segment=[1 segment];
[segu,segui]=unique(segment);

for k=1:length(segu)
    iu=find(segment==segu(k));
    imin=min(iu);
    imax=max(iu);
    
    bef=max(iu(1)-1,1);
    
    slope=(round((conduc(iu(end))-conduc(iu(1)))*1000)/1000)./(S_FILE.DATES(iu(end))-S_FILE.DATES(iu(1)));
    disp('------------------------------------------------------------------')
    disp('OW correction')
    disp('-----------------------------')
    texte=['cycles ' num2str(prof_no(imin)) ' to ' num2str(prof_no(imax)) ': Correction Y(' num2str(prof_no(imin)) ') ' num2str((conduc(bef))) ' +/- ' num2str((conduc_err(bef))) ' and Slope ' num2str(slope) '/yr' ];
    disp(texte)
    disp('------------------------------------------------------------------')
end
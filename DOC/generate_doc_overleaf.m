% -========================================================
%   USAGE : generate_doc(tabfloat,tabdac,configow,varargin)
%   PURPOSE : run OWC software for a float or a group of floats, using a specified OW config
% -----------------------------------
%   INPUT :
%    tabfloat  (char or cell of chars -size n_floatsx1)    e.g. '6900258' or {'6900258', '3901954'}
%    tabdac    (char or cell of chars -size n_floatsx1)    e.g. 'coriolis' or {'coriolis', 'bodc'}
%    configow  (floats or cell of floats -size n_floatsx1) e.g.  149       or {149,149}      % config number ow
%
%   OPTIONNAL INPUT :
%    PROFREF   (array or cell of array -size n_floatsx1)  e.g. 1 or {1,[],10}  float profiles for which reference profile is plotted
%    PROFREFADJ(array or cell of array -size n_floatsx1)  e.g. 1 or {1,[],10}  float adjusted profiles for which reference profile  is plotted
%    PROFFLAG  (array or cell of array -size n_floatsx1)  e.g. 1 or {1 or '1','1D','2'}  float  profiles for  which flag are plotted
%    PROFCLOSE (array or cell of array -size n_floatsx1)  e.g. 1 or {1,[],10}  float profiles for which closest Argo profile is plotted
%    COMP_GDAC  (logical) 1 if comparison to GDAC profile is plotted , 0 if not (default)
%    PLOTDEPTHDEP (logical) 1 if plot of depth dependance is included (0 if not, default)
%    TITLE     (char) title of the document  (default: 'Delayed mode analysis of salinity data acquired by Argo floats';)
%    SUBTITLE  (char)  subtilte of the document (default: float_name)  e.g.: 'deployed in 2017 (PI : G. MAZE)'  e.g '6902757 (PI : G. MAZE)' )
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created (2016) ccabanes
%
%   CALLED SUBROUTINES:
% -------------------------------------
%
% ========================================================
function generate_doc_overleaf(tabfloat,tabdac,configow,varargin)

%init_path

if iscell(tabfloat)==0;tabfloat=cellstr(tabfloat);end
if iscell(tabdac)==0;tabdac=cellstr(tabdac);end
if iscell(configow)==0;configow={configow};end
if length(tabfloat)>1&length(tabdac)==1
    tabdac=repmat(tabdac,1,length(tabfloat));
end
if length(tabfloat)>1&length(configow)==1
    configow=repmat(configow,1,length(tabfloat));
end


float_list=tabfloat;
for k=1:length(tabfloat)
    Num_Config{k}=num2str(configow{k});
end

CONF = load_configuration('config.txt');
DIR_FTP=CONF.DIR_FTP;
DIR_FTP_ORIG=CONF.DIR_FTP_CORIOLIS;
DIR_PLOT=CONF.DIR_PLOT;
DIR_DATA=CONF.DIR_DATA;
DIR_PLOTOW=[CONF.DIR_DATA 'float_plots/CONFIG'];
DIR_DMQC=CONF.DIR_DM_FILES;
DIR_CODES=CONF.DIR_CODES;
AUTHORS=CONF.REPORT_AUTHORS;

% INPUT PARAMETERS
n=length(varargin);
if n/2~=floor(n/2)
    error('check the imput arguments')
end
f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);
% default
Profrefadj=repmat({[]},1,length(tabfloat));
Profref=repmat({[]},1,length(tabfloat));
Profclose=repmat({[]},1,length(tabfloat));
Profflag=repmat({[]},1,length(tabfloat));
TITLE='Delayed mode analysis of salinity data acquired by Argo floats ';

if length(unique(tabfloat))>1 ;
    SUBTITLE=['Floats ' strjoin(unique(tabfloat), ', ')];
else
    M = read_netcdf_allthefile([DIR_FTP tabdac{1} '/' float_list{1} '/' float_list{1} '_meta.nc']);
    pi_name=deblank(M.pi_name.data');
    pi_name_red=reduce_pi_name(pi_name);
    SUBTITLE=['Float ' strjoin(unique(tabfloat) , ', ') ' (' pi_name_red ')' ];
end

COMP_GDAC=0;
PLOTDEPTHDEP=0;
FORCE_DEEP=0;
NEW_CPCOR=-11.6e-8;

if isfield(s,'PROFREF')==1;Profref=s.PROFREF;end;
if isfield(s,'PROFREFADJ')==1;Profrefadj=s.PROFREFADJ;end;
if isfield(s,'PROFFLAG')==1;Profflag=s.PROFFLAG;end;
if isfield(s,'PROFCLOSE')==1;Profclose=s.PROFCLOSE;end
if isfield(s,'TITLE')==1;TITLE=s.TITLE;end;
if isfield(s,'SUBTITLE')==1;SUBTITLE=s.SUBTITLE;end;
if isfield(s,'COMP_GDAC')==1;COMP_GDAC=s.COMP_GDAC;end;
if isfield(s,'PLOTDEPTHDEP')==1;PLOTDEPTHDEP=s.PLOTDEPTHDEP;end;
if isfield(s,'FORCE_DEEP')==1;FORCE_DEEP=s.FORCE_DEEP;end;
if isfield(s,'NEW_CPCOR')==1;NEW_CPCOR=s.NEW_CPCOR;end;

if iscell(Profref)==0;Profref={Profref};end
if iscell(Profrefadj)==0;Profrefadj={Profrefadj};end
if iscell(Profclose)==0;Profclose={Profclose};end
if iscell(Profflag)==0;Profflag={Profflag};end


if length(unique(Num_Config))>1 ;
    str_conf=['configurations ' strjoin(unique(Num_Config), ', ')];
else
    str_conf=['configuration ' strjoin(unique(Num_Config), ', ')];
end


for ik=1:length(float_list)
    M = read_netcdf_allthefile([DIR_FTP tabdac{ik} '/' float_list{ik} '/' float_list{ik} '_meta.nc']);
    [isfound]=findstr_tab(M.config_parameter_name.data,'ParkPressure_dbar');
    if sum(isfound)==1
        park_press=median(M.config_parameter_value.data(:,isfound));park_press_def=0;
        
    else
        park_press=1000;park_press_def=1;
    end
    [isfound]=findstr_tab(M.config_parameter_name.data,'ProfilePressure_dbar');
    if sum(isfound)==1
        prof_press=median(M.config_parameter_value.data(:,isfound));prof_press_def=0;
    else
        prof_press=2000;prof_press_def=1;
    end
end


% Added By CC+TR 25.09.20
dir_tex='OVERLEAF/';
if exist('OVERLEAF','dir')
    [success,message,messageid]=rmdir('OVERLEAF','s');
    if success==0
        warning('OVERLEAF DIRECTORY CANNOT BE REMOVED')
    else
        mkdir(dir_tex);
    end
else
    mkdir(dir_tex);
end

% lit le fichier des corrections
num=[];

file_cor='correction_float.csv';

if exist(file_cor)==0
    file_cor='correction_float_template.csv';
    [num,wmo_corr,col1,col2_1,col2_2,col3,col4,col5]=get_txtfile_col(file_cor,';');
    wmo_corr=float_list;
else
    
    [num,wmo_corr,col1,col2_1,col2_2,col3,col4,col5]=get_txtfile_col(file_cor,';');
    %copyfile(file_cor ,[dir_tex file_cor]);
end


fw1=fopen(['./OVERLEAF/ebauche_rapport.tex'],'w');

% ecriture de l'entete du document + titre
%-----------------------------------------
file_cor
write_titletex(fw1,TITLE, SUBTITLE, AUTHORS,file_cor,float_list );

fprintf(fw1,'%s\n', [' ']);
fprintf(fw1,'%s\n', [' ']);
fprintf(fw1,'%s\n', ['%------------------------------------------------------ ']);
fprintf(fw1,'%s\n', [' ']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fw1,'%s\n', ['\section{Presentation}']);
fprintf(fw1,'%s\n', [' ']);

if length(tabfloat)>1
    str=['Delayed Mode analysis was performed for each  float (see table \ref{tab1}). First, salinity and temperature profiles were visually checked and compared to nearby reference profiles using verif\_flag programs when necessary. Real time QC flags were verified and modified if necessary (see table \ref{tab2}). The OWC method was then run to estimate a salinity offset or/and a salinity drift, using, if possible, historical  CTD or Argo profiles as reference databases. Finally, corrections were applied in the netcdf files when we though it was necessary(see table \ref{tab4}).'];
else
    str=['Delayed Mode analysis was performed float ' tabfloat{1} '. First, salinity and temperature profiles were visually checked and compared to nearby reference profiles using verif\_flag programs when necessary. Real time QC flags were verified and modified if necessary (see table \ref{tab2}). The OWC method was then run to estimate a salinity offset or/and a salinity drift, using, if possible, historical CTD or Argo profiles as reference databases. Finally, corrections were applied in the netcdf files when we though it was necessary(see table \ref{tab4}).'];
end

fprintf(fw1,'%s\n', str);
%%----------------------------------------------------------------------------
% ecriture de la partie "table 1" WMo -Launch date  -Centre PI - Last cycle analysed
%----------------------------------------------------------------------------
% preparation du tableau
fprintf(fw1,'%s\n', ['\begin{table}[h]']);
fprintf(fw1,'%s\n', ['$$']);
%fprintf(fw1,'%s\n', ['\begin{tabular}{|l|c|c|c|c|c|}']);
fprintf(fw1,'%s\n', ['\begin{tabular}{|l|c|c|c|c|}']);
fprintf(fw1,'%s\n', ['\hline']);
%fprintf(fw1,'%s\n', ['WMO Number & Launch date & Centre & PI & Last cycle analysed &Cycle Duration \\']);
fprintf(fw1,'%s\n', ['WMO Number & Launch date & Centre & PI & Last cycle analysed  \\']);

%fprintf(fw1,'%s\n', ['& & & & (Active/NotActive) &  \\']);
fprintf(fw1,'%s\n', ['& & & & (Active/NotActive)   \\']);

fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['\hline']);
%--------------------

% remplissage du tableau
NcVar.juld.name='JULD';
NcVar.config_mission_number.name='CONFIG_MISSION_NUMBER';
NcVar.cycle_number.name='CYCLE_NUMBER';
NcVar.date_update.name='DATE_UPDATE';
disp('Recherche les info dans les fichiers meta et prof')

for ik=1:length(float_list)
    % Table 1 recherche les info dans les fichiers meta et prof
    M = read_netcdf_allthefile([DIR_FTP tabdac{ik} '/' float_list{ik} '/' float_list{ik} '_meta.nc']);
    P = read_netcdf_allthefile([DIR_FTP tabdac{ik} '/' float_list{ik} '/' float_list{ik} '_prof.nc'],NcVar);
    dir_upload=dir([DIR_FTP_ORIG tabdac{ik} '/' float_list{ik} '/'] ); 
    upload_date=dir_upload(find(findstr_tab({dir_upload.name},'profiles'))).datenum; % date du dernier load. On verifie si le flotteur est actif ou pas a cette date
    if isempty(strfind(DIR_FTP_ORIG, '/home/coriolis_exp/spool/co05/co0508/'))==0; % local copy of GDAC fts server
        upload_date=datenum(date);
    end
    pi_name=deblank(M.pi_name.data');
    
    pi_name_red=reduce_pi_name(pi_name);
    
    thedate=datestr(datevec(M.launch_date.data','YYYYmmddHHMMSS'),'dd/mm/YYYY');
    ikl=find(findstr_tab(M.config_parameter_name.data,'CONFIG_CycleTime_days'));
    
    if isempty (ikl)==0
        duration=(M.config_parameter_value.data(:,ikl));
        
    else
        ikl=find(findstr_tab(M.config_parameter_name.data,'CONFIG_CycleTime_hours'));
        if isempty(ikl)==0
            duration=(M.config_parameter_value.data(:,ikl))/24;
        end
        
    end
    if (P.juld.data(end)+20)<( upload_date-datenum('19500101','YYYYmmdd'))
        str_act='NA';
    else
        str_act='A';
    end
    
    if isempty(ikl)==0
        P.duration=duration(P.config_mission_number.data);
        if (P.juld.data(end)+3*duration(end))<( upload_date-datenum('19500101','YYYYmmdd'))
            str_act='NA';
        else
            str_act='A';
        end
        
%         str_duration=[ 'cy.' num2str(P.cycle_number.data(1)) '-'];
%         ico=1;
        %for lm=1:length(P.duration)-1
            %if P.duration(lm)~=P.duration(lm+1)
                
                %str_duration=[str_duration  num2str(P.cycle_number.data(lm)) ': ' num2str(P.duration(lm)) ' days'];
               % if ico<2;
                    %fprintf(fw1,'%s\n', [float_list{ik} ' & ' thedate '&' M.data_centre.data' '&' pi_name_red '&' num2str(P.cycle_number.data(end)) '(' str_act ') &' str_duration  '\\']);
                    fprintf(fw1,'%s\n', [float_list{ik} ' & ' thedate '&' M.data_centre.data' '&' pi_name_red '&' num2str(P.cycle_number.data(end)) '(' str_act ') \\']);
               % else
                    %fprintf(fw1,'%s\n', [' & & & & &' str_duration  '\\']);
                %end
%                 str_duration=['cy.' num2str(P.cycle_number.data(lm+1)) '-'];
%                 ico=ico+1;
           % end
            
        %end
%         str_duration=[str_duration  num2str(P.cycle_number.data(lm+1)) ': ' num2str(P.duration(lm+1)) ' days'];
%     else
%         str_duration='-';
%         ico=1;
    end
%     if ico<2;
%         fprintf(fw1,'%s\n', [float_list{ik} ' & ' thedate '&' M.data_centre.data' '&' pi_name_red '&' num2str(P.cycle_number.data(end)) '(' str_act ') &' str_duration  '\\']);
%     else
%         fprintf(fw1,'%s\n', [' & & & & &' str_duration  '\\']);
%     end
    fprintf(fw1,'%s\n', ['\hline']);
    
   
    
    
end

% fin du tableau et legende

%fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['\end{tabular}']);
fprintf(fw1,'%s\n', ['$$']);
fprintf(fw1,'%s\n', ['\caption{Information on the floats analysed}']);
fprintf(fw1,'%s\n', ['\label{tab1}']);
fprintf(fw1,'%s\n', ['\end{table}'])
% fin de la table 1
%----------------------------------------------------------------------------
fprintf(fw1,'%s\n', ['%------------------------------------------------------ ']);
fprintf(fw1,'%s\n', [' ']);



%%----------------------------------------------------------------------------
% ecriture d'un tableau mission
%----------------------------------------------------------------------------
% preparation du tableau
fprintf(fw1,'%s\n', ['\begin{table}[h]']);
fprintf(fw1,'%s\n', ['$$']);
fprintf(fw1,'%s\n', ['\begin{tabular}{|l|c|c|c|c|c|}']);
fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['WMO Number & Cycles & Cycle Duration & Park Pressure & Profile Pressure & ISA Ice Detection  \\']);
fprintf(fw1,'%s\n', ['&  & (days) & (db) & (db) & (degC) \\']);
fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['\hline']);

for ik=1:length(float_list)
    Mf=read_file_meta([DIR_FTP tabdac{ik} '/' float_list{ik} '/' float_list{ik} '_meta.nc']);
    P = read_netcdf_allthefile([DIR_FTP tabdac{ik} '/' float_list{ik} '/' float_list{ik} '_prof.nc'],NcVar);
    ico=1;
    str_cy=[ 'cy.' num2str(P.cycle_number.data(1)) '-'];
    for lm=1:length(P.config_mission_number.data)-1
        if P.config_mission_number.data(lm)~=P.config_mission_number.data(lm+1)
            str_cy=[ str_cy num2str(P.cycle_number.data(lm))  ];
            mission=P.config_mission_number.data(lm);
            fmis=find(Mf.config_mission_number==mission);
            if length(fmis)==1;
                str_day=num2str(Mf.CycleTime(fmis)); str_day=strrep(str_day,'NaN','-');
                str_park=num2str(Mf.ParkPressure(fmis));str_park=strrep(str_park,'NaN','-');
                str_prof=num2str(Mf.ProfilePressure(fmis));str_prof=strrep(str_prof,'NaN','-');
                str_ice=num2str(Mf.IceDetection(fmis));str_ice=strrep(str_ice,'NaN','-');
            else
                str_park='-';
                str_prof='-';
                str_ice='-';
            end
            if ico<2;
                fprintf(fw1,'%s\n', [float_list{ik} ' &' str_cy ' & ' str_day ' & ' str_park ' & ' str_prof ' & ' str_ice ' \\']);
            else
                fprintf(fw1,'%s\n', [' & ' str_cy ' & ' str_day ' & ' str_park ' & ' str_prof ' & ' str_ice ' \\']);
            end
            str_cy=['cy.' num2str(P.cycle_number.data(lm+1)) '-'];
            ico=ico+1;
        end
        
    end
    
    str_cy=[str_cy  num2str(P.cycle_number.data(lm+1)) ];
    mission=P.config_mission_number.data(lm+1);
    fmis=find(Mf.config_mission_number==mission);
    if length(fmis)==1;
        str_day=num2str(Mf.CycleTime(fmis));str_day=strrep(str_day,'NaN','-');
        str_park=num2str(Mf.ParkPressure(fmis));str_park=strrep(str_park,'NaN','-');
        str_prof=num2str(Mf.ProfilePressure(fmis));str_prof=strrep(str_prof,'NaN','-');
        str_ice=num2str(Mf.IceDetection(fmis));str_ice=strrep(str_ice,'NaN','-');
    else
        str_park='-';
        str_prof='-';
        str_ice='-';
    end
    if ico<2;
       fprintf(fw1,'%s\n', [float_list{ik} ' &' str_cy ' & ' str_day ' & ' str_park ' & ' str_prof ' & ' str_ice ' \\']);
    else
        fprintf(fw1,'%s\n', [' & ' str_cy ' & ' str_day ' & ' str_park ' & ' str_prof ' & ' str_ice ' \\']);
    end
    fprintf(fw1,'%s\n', ['\hline']);
end
%fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['\end{tabular}']);
fprintf(fw1,'%s\n', ['$$']);
fprintf(fw1,'%s\n', ['\caption{Configuration Parameters for each mission }']);
fprintf(fw1,'%s\n', ['\label{tab1b}']);
fprintf(fw1,'%s\n', ['\end{table}'])
% fin de la table 1
%----------------------------------------------------------------------------
fprintf(fw1,'%s\n', ['%------------------------------------------------------ ']);
fprintf(fw1,'%s\n', [' ']);
fprintf(fw1,'%s\n', '\clearpage');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% DMQC SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fw1,'%s\n', ['\section{DMQC Summary}']);

fprintf(fw1,'%s\n', ['\subsection{Verification of RT QC flags}']);


str=['Real Time QC flags were verified and modified if necessary. Table \ref{tab2} gives the list of flags that have been modified during the delayed mode process.'];

fprintf(fw1,'%s\n', str);

%----------------------------------------------------------------------------
%  ecriture de la partie TABLE2  : Correction des flags
%----------------------------------------------------------------------------
% préparation du tableau
fprintf(fw1,'%s\n', ['\setlongtables']);
fprintf(fw1,'%s\n', ['\begin{longtable}{|l|c|c|c|c|c|c|}']);
%fprintf(fw1,'%s\n', ['\begin{table}[h!]']);
%fprintf(fw1,'%s\n', ['$$']);
%fprintf(fw1,'%s\n', ['\begin{tabular}{|l|c|c|c|c|c|c|}']);
fprintf(fw1,'%s\n', ['\nobreakhline']);
fprintf(fw1,'%s\n', ['WMO Number & Cycle & Param & Old flag & New flag & Levels & Date of modification \\']);
fprintf(fw1,'%s\n', ['\nobreakhline']);
fprintf(fw1,'%s\n', ['\nobreakhline']);
fprintf(fw1,'%s\n', ['\endhead']);
fprintf(fw1,'%s\n', ['\endfoot']);

NcVar.history_qctest.name=upper('history_qctest');
NcVar.history_date.name=upper('history_date');
NcVar.history_software.name=upper('history_software');
NcVar.history_software_release.name=upper('history_software_release');
NcVar.history_parameter.name=upper('history_parameter');
NcVar.history_previous_value.name=upper('history_previous_value');
NcVar.history_start_pres.name=upper('history_start_pres');
NcVar.history_stop_pres.name=upper('history_stop_pres');
NcVar.platform_number.name=upper('platform_number');
NcVar.cycle_number.name=upper('cycle_number');
NcVar.history_action.name=upper('history_action');
NcVar.direction.name='DIRECTION';
NcVar.temp_qc.name='TEMP_QC';
NcVar.psal_qc.name='PSAL_QC';
NcVar.pres.name='PRES';

disp('Table 2: recherche des informations dans l''history')
% Recherche des informations dans l'history
for ik=1:length(float_list)
    %[Co,Dim]=create_multi_from_mono(DIR_FTP,float_list{ik},tabdac{ik},'CR','Primary sampling',NcVar);
    [file_list]=select_float_files_on_ftp(float_list{ik},tabdac{ik},DIR_FTP,'C');
    [Co,Dim]=create_multi_from_filelist(float_list{ik},tabdac{ik},DIR_FTP,file_list,'Primary sampling',NcVar);
    %(floatname,dacname,DIRFTP,file_list,vertical_sampling_scheme,Param);
    %keyboard
    %[Co,Dim]=create_multi_from_mono(DIR_DMQC,float_list{ik},'','C','Primary Sampling',NcVar);
    read_history_cond(Co,Dim,NcVar,fw1) % lit l'historique et ecrit dans la table
    Co=replace_fill_bynan(Co);
    mincy=1;
    maxcy=max(Co.cycle_number.data);
    cy_missing{ik}=setdiff([mincy:maxcy],Co.cycle_number.data');
    
    % recupere quelques indications sur les niveaux de pression:
    minpres=min(Co.pres.data');
    maxpres=max(Co.pres.data');
    isnotdive=find(maxpres<50);
    cy_isnotdive{ik}=Co.cycle_number.data(isnotdive);
    level_lack=find(minpres>20);
    cy_level_lack{ik}=Co.cycle_number.data(level_lack);
    cy_level_lack_value{ik}=minpres(level_lack);
    tab(ik,1)=Co.cycle_number.data(1);
    tab(ik,2)=Co.cycle_number.data(end);
end

fprintf(fw1,'%s\n', ['\caption{Modified flags during DM analysis}']);
fprintf(fw1,'%s\n', ['\label{tab2}']);
fprintf(fw1,'%s\n', ['\end{longtable}']);
% fin de la table 2
%----------------------------------------------------------------------------


fprintf(fw1,'%s\n', ['%------------------------------------------------------ ']);
fprintf(fw1,'%s\n', [' ']);

% ----------------------------------------------------------------------------
% note les cycles ou des inversions de densité sont trouvees et differents problèmes (niveaux manquants, flotteur qui n'a pas plongé)
% ----------------------------------------------------------------------------
str=['For each float, we report here the list of cycles for which a density inversion was detected in real time (with a treshold value of 0.03). This sometimes reveals a problem with the conductivity sensor and it is necessary to particularly check these profiles in delayed time. Moreover, when density inversion are flagged in RT, it is often necessary to modified flags in DM: often, the temperature does not need to be flagged at 4 and  not all the salinity measurements flagged in RT need a flag 4. We also report here some anomalies e.g. a float that did not dive for a given cycle or missing cycles. '];

fprintf(fw1,'%s\n', str);

fprintf(fw1,'%s\n', ['\begin{itemize}']);

disp('Recherche des inversions de densite dans les logs')
comptitem=0;
% lit les log des inversions de densite
for ik=1:length(float_list)
    file_inv_name=[DIR_PLOT 'density_anomaly/' float_list{ik} '/' float_list{ik} '_chkinv_sigflagnotused.txt'];
    if exist(file_inv_name)
        fwr=fopen(file_inv_name,'r');
        C=textscan(fwr,'%s','Delimiter','\n');
        fclose(fwr)
        kinv=0;
        str_to_write='';
        s1='';
        clear strcy
        for ll=1:length(C{1})
            if findstr_tab(C{1}{ll},'Cycle')==1;
                kinv=kinv+1;
                strcy{kinv}=(strtok(C{1}{ll},'Cycle '));
            end
        end
        
        if kinv>0;
            s1=strjoin(unique(strcy,'stable'),', ');
        end
        if isempty(s1)==0;
            if length(unique(strcy))>1 ;
                strcycle='cycles';
            else
                strcycle='cycle';
            end
            str_to_write=[ str_to_write float_list{ik} ' - Density inversions are found ' strcycle ': ' s1 '. '];
        else
            str_to_write=[ str_to_write float_list{ik} ' - No Density inversions. '];
        end
        if isempty(cy_isnotdive{ik})==0;
            if length(cy_isnotdive{ik})>1;
                strcycle='cycles';
            else
                strcycle='cycle';
            end
            str_to_write=[ str_to_write 'The float did not dive at ' strcycle ': ' strjoin(cellstr(num2str(cy_isnotdive{ik}))',', ') '. '];
        end
        if isempty(cy_level_lack{ik})==0; fprintf(fw1,'%s\n', '');
            if length(cy_level_lack{ik})>1;
                strcycle='cycles ';
            else
                strcycle='cycle ';
            end
            str_to_write=[ str_to_write 'Upper levels date are missing ' strcycle ':' strjoin(cellstr(num2str(cy_level_lack{ik}'))',', ') '. '];
        end
        if isempty(cy_missing{ik})==0; fprintf(fw1,'%s\n', '');
            if length(cy_missing{ik})>1;
                strcycle='cycles ';
            else
                strcycle='cycle ';
            end
            str_to_write=[ str_to_write 'Missing ' strcycle ':' strjoin(cellstr(num2str(cy_missing{ik}'))',', ') '. '];
        end
        if isempty(str_to_write)==0;
            comptitem=comptitem+1;
            fprintf(fw1,'%s\n', ['\item ' str_to_write]);
        end
    end
end
if comptitem==0
    %fprintf(fw1,'%s\n', ['\item No density inversion was detected for the floats']);
end
fprintf(fw1,'%s\n',['\end{itemize}']);

fprintf(fw1,'%s\n', '\clearpage');

%% ----------------------------------------------------------------------------
%     TABLE 3 : Salinité correction applied
%%----------------------------------------------------------------------------

fprintf(fw1,'%s\n', ['\subsection{Salinity corrections applied}']);


dir_fig=[ DIR_PLOT 'verif_profil1/' float_list{ik} '/'];% Added By TR 28.09.20
cpcor_filename=['CPCOR_analysis_' float_list{ik}   '_1.png'];
if exist([dir_fig,cpcor_filename])
    copyfile([dir_fig,cpcor_filename],[dir_tex,cpcor_filename]);
end

if exist([dir_tex,cpcor_filename])==0&FORCE_DEEP==0
    % lit le fichier texte des corrections appliquees,
    if exist(file_cor)~=0
        %[num,wmo_corr,col1,col2_1,col2_2,col3,col4,col5]=get_txtfile_col(file_cor,';');
        
        if sum(ismember(wmo_corr,float_list))>=1 % on fait un tableau
            fprintf(fw1,'%s\n', ['\renewcommand\arraystretch{1.2}']);
            fprintf(fw1,'%s\n', ['\begin{table}[h]']);
            fprintf(fw1,'%s\n', ['$$']);
            fprintf(fw1,'%s\n', ['\begin{tabular}{|l|m{4cm}|m{6cm}|m{5cm}|}']);
            fprintf(fw1,'%s\n', ['\hline']);
            fprintf(fw1,'%s\n', ['           &  \multicolumn{2}{c|}{Calibration}    &      \\']);
            
            fprintf(fw1,'%s\n', ['   WMO        & \center{Comparison with the reference CTD cast} & \center{Correction from OWC   method}                &   Correction applied in the D files       \\']);
            fprintf(fw1,'%s\n', ['  Number      &            &          &                         \\']);
            fprintf(fw1,'%s\n', ['\hline']);
            fprintf(fw1,'%s\n', ['\hline']);
            for ik=1:length(float_list)
                iil=find(findstr_tab(wmo_corr,float_list{ik}));
                if isempty(iil)==0
                    fprintf(fw1,'%s\n', [float_list{ik} '        &   ' col1{iil}  '                 & ' [col2_1{iil} ' ' col2_2{iil}  ' (config. ' Num_Config{ik}   ')' ]     '            & ' col3{iil}        '          \\']);
                    %fprintf(fw1,'%s\n', ['               &                     & '    col2_2{iil}  ' (config. ' Num_Config{ik}   ')'                          '  &                           \\']);
                    fprintf(fw1,'%s\n', '\hline');
                    fprintf(fw1,'%s\n', '%------------------------------------------------------------------------------------------------------------');
                end
            end
            
            fprintf(fw1,'%s\n', '\end{tabular}');
            fprintf(fw1,'%s\n', '$$');
            fprintf(fw1,'%s\n', '\caption{Salinity corrections for the  floats proposed by the OWC method or by comparison with a shipboard CTD reference profile. Uncertainties are the statistical uncertainties from the OWC method.}');
            fprintf(fw1,'%s\n', '\label{tab4}');
            fprintf(fw1,'%s\n', '\end{table}');
        end
    end
else % si le CPCOR a ete corrige
    % lit le fichier texte des corrections appliquees,
    if exist(file_cor)~=0
        %[num,wmo_corr,col1,col2_1,col2_2,col3,col4,col5]=get_txtfile_col(file_cor,';');
        
        if sum(ismember(wmo_corr,float_list))>=1 % on fait un tableau
            fprintf(fw1,'%s\n', ['\renewcommand\arraystretch{1.2}']);
            fprintf(fw1,'%s\n', ['\begin{table}[h]']);
            fprintf(fw1,'%s\n', ['$$']);
            fprintf(fw1,'%s\n', ['\begin{tabular}{|l|m{2cm}|m{3cm}|m{6cm}|m{4cm}|}']);
            fprintf(fw1,'%s\n', ['\hline']);
            fprintf(fw1,'%s\n', ['      &     &  \multicolumn{3}{c|}{Calibration (with new CPcorr value applied)}         \\']);
            
            fprintf(fw1,'%s\n', ['   WMO  &  new CPcorr     & \center{Comparison with the reference CTD cast} & \center{Correction from OWC   method}                &   Correction applied in the D files       \\']);
            fprintf(fw1,'%s\n', ['  Number    &              &            &          &                         \\']);
            fprintf(fw1,'%s\n', ['\hline']);
            fprintf(fw1,'%s\n', ['\hline']);
            for ik=1:length(float_list)
                iil=find(findstr_tab(wmo_corr,float_list{ik}));
                if isempty(iil)==0
                    fprintf(fw1,'%s\n', [float_list{ik} '        & $' num2str(NEW_CPCOR) '$ & ' col1{iil}  '                 & ' [col2_1{iil} ' ' col2_2{iil}  ' (config. ' Num_Config{ik}   ')' ]     '            & ' col3{iil}        '          \\']);
                    %fprintf(fw1,'%s\n', ['               &                     & '    col2_2{iil}  ' (config. ' Num_Config{ik}   ')'                          '  &                           \\']);
                    fprintf(fw1,'%s\n', '\hline');
                    fprintf(fw1,'%s\n', '%------------------------------------------------------------------------------------------------------------');
                end
            end
            
            fprintf(fw1,'%s\n', '\end{tabular}');
            fprintf(fw1,'%s\n', '$$');
            fprintf(fw1,'%s\n', '\caption{Salinity corrections for the  floats proposed by the OWC method or by comparison with a shipboard CTD reference profile once the new Cpcorr  value has been applied to the conductivity data. Uncertainties are the statistical uncertainties from the OWC method.}');
            fprintf(fw1,'%s\n', '\label{tab4}');
            fprintf(fw1,'%s\n', '\end{table}');
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES POUR CHAQUE FLOTTEUR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Début des figures
for ik=1:length(float_list)
    labell=['_' num2str(ik)];
    fprintf(fw1,'%s\n', ' ');
    fprintf(fw1,'%s\n', ' ');
    fprintf(fw1,'%s\n', ['%------------------------------------------------------ ']);
    
    fprintf(fw1,'%s\n', '\clearpage');
    
    %% figure des déplacements
    fprintf(fw1,'%s\n', ['\section {Float ' float_list{ik} '}']);
    fprintf(fw1,'%s\n', ['\subsection {Trajectory}']);
    
    fprintf(fw1,'%s\n', ' %%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf(fw1,'%s\n', '\begin{figure}[h!]');
    fprintf(fw1,'%s\n', '\begin{flushleft}');
   
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.28\linewidth}');
    %fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/'];% Added By TR 25.09.20
    file =[float_list{ik} '_pos2_flagnotused.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    str=['\includegraphics[width=5cm,trim= 10 10 10 10, clip=true]{' file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{1\linewidth}');
    %fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/'];% Added By TR 25.09.20
    file =[float_list{ik} '_pos_flagnotused.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %file =[ DIR_PLOT 'preliminaire/' float_list{ik} '/' float_list{ik} '_pos_flagnotused.png'];
    str=['\includegraphics[width=15cm,trim= 10 10 10 10, clip=true]{' file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    %   fprintf(fw1,'%s\n', '$$');
    str=['Float ' float_list{ik} '. Trajectory of the float and bathymetry. Parking depth is: ' num2str(park_press) 'm and profile depth is: ' num2str(prof_press) 'm. Bathymetric contours at  float''s parking depth $\pm$ 30m are plotted in green, bathymetric contours at  float''s profile depth $\pm$ 30m are plotted in red, bathymetric contours between profile depth and parking depth are plotted every 200m in magenta and bathymetric contours between parking depth and surface are plotted every 200m in blue. )'];
   % str=['Float ' float_list{ik} '. Trajectory of the float and bathymetry. Parking depth is: ' num2str(park_press) 'm and profile depth is: ' num2str(prof_press) 'm. Bathymetric contours at  float''s parking depth $\pm$ 30m are plotted in green, bathymetric contours at  float''s profile depth $\pm$ 30m are plotted in red, bathymetric contours between profile depth and parking depth are plotted every 200m in magenta and bathymetric contours between parking depth and surface are plotted every 200m in blue. Grey arrows are for POSITION\_QC =1,2 (good, probably good); magenta arrows are for POSITION\_QC=3,4 (probably bad, bad); white arrows are for POSITION\_QC=8 (interpolated)'];

    fprintf(fw1,'%s\n', ['\caption{' str '}']);
    fprintf(fw1,'%s\n', ['\label{fig1' labell '}']);
    fprintf(fw1,'%s\n', '\end{flushleft}');
    fprintf(fw1,'%s\n', '\end{figure}');
    fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf(fw1,'%s\n', ' ');
    
    fprintf(fw1,'%s\n', '\clearpage');
    fprintf(fw1,'%s\n', ['\subsection {Sections along the float trajectory - raw data}']);
    
    %% figure des sections TPOT,PSAL,SIG0
    
    fprintf(fw1,'%s\n', ' %%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf(fw1,'%s\n','\begin{figure}[h!]');
    fprintf(fw1,'%s\n','\begin{center}');
    
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.48\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/' ];% Added By TR 25.09.20
    file =[float_list{ik} '_TPOT_interp_flagnotused.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{' file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.48\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/'];% Added By TR 25.09.20
    file =[float_list{ik} '_SIG0_interp_flagnotused.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=270,origin=c,width=\textwidth,trim= 15 10 20 45, clip=true]{' DIR_PLOT 'preliminaire/' float_list{ik} '/' float_list{ik} '_SIG0_interp_flagnotused.png}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{' file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.8\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/'];% Added By TR 25.09.20
    file =[float_list{ik} '_PSAL_interp_flagnotused.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=270,origin=c,width=\textwidth,trim= 15 10 20 45, clip=true]{' DIR_PLOT 'preliminaire/' float_list{ik} '/' float_list{ik} '_PSAL_interp_flagnotused.png}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{' file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    str=['Float ' float_list{ik} '. Potential temperature, Sig0 and salinity sections along the float trajectory (raw data, flags not used)'];
    fprintf(fw1,'%s\n', ['\caption{' str '}']);
    fprintf(fw1,'%s\n', ['\label{fig2' labell '}']);
    fprintf(fw1,'%s\n', '\end{center}');
    fprintf(fw1,'%s\n', '\end{figure}');
    fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf(fw1,'%s\n', ' ');
    
     %% figure des theta/S
     
    fprintf(fw1,'%s\n', '\clearpage');
    fprintf(fw1,'%s\n', ['\subsection {Theta/S diagrams - raw data}']);
    % figure du diag OS
    fprintf(fw1,'%s\n', ' %%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf(fw1,'%s\n','\begin{figure}[h!]');
    fprintf(fw1,'%s\n','\begin{center}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.65\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/'];% Added By TR 25.09.20
    file =[float_list{ik} '_thetaS_flagnotused.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=270,origin=c,width=\textwidth,trim= 15 10 20 50, clip=true]{' DIR_PLOT 'preliminaire/' float_list{ik} '/' float_list{ik} '_thetaS_flagnotused.png}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{' file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.65\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/'];% Added By TR 25.09.20
    file =[float_list{ik} '_thetaS_zoom_flagnotused.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=270,origin=c,width=\textwidth,trim= 15 10 20 50, clip=true]{' DIR_PLOT 'preliminaire/' float_list{ik} '/' float_list{ik} '_thetaS_zoom_flagnotused.png}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{' file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    str=['Float ' float_list{ik} '. Theta/S diagrams of the raw data, with the potential temperature referenced to 0db. Full profiles (upper panel) and zoom below 1500m (lower panel). Flags are not used'];
    fprintf(fw1,'%s\n', ['\caption{' str '}']);
    fprintf(fw1,'%s\n', '\label{fig3}');
    fprintf(fw1,'%s\n', '\end{center}');
    fprintf(fw1,'%s\n', '\end{figure}');
    
    %% Figure des données techniques
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/'];% Added By TR 28.09.20
    filenamecmp=[float_list{ik} '_surface_pres.pdf'];
    if exist([dir_fig,filenamecmp])
        copyfile([dir_fig,filenamecmp],[dir_tex,filenamecmp]);
    end
    
    if exist([dir_tex,filenamecmp])
        fprintf(fw1,'%s\n', '\clearpage');
        fprintf(fw1,'%s\n', ['\subsection {Technical data : surface pressure - battery - pump or valve actions}']);
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\begin{figure}[h!]');
        fprintf(fw1,'%s\n', '\begin{center}');
        fprintf(fw1,'%s\n', '$$');
        str=['\includegraphics[width=15cm,trim=  20 80 20 80, angle=0,clip=true]{'  filenamecmp '}'];
        fprintf(fw1,'%s\n', str);
        fprintf(fw1,'%s\n', '$$');
        str=['Float ' float_list{ik} ': Some technical data as read in the technical file'];
        fprintf(fw1,'%s\n', ['\caption{' str '}']);
        fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_flag}']);
        fprintf(fw1,'%s\n', '\end{center}');
        fprintf(fw1,'%s\n', '\end{figure}');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        
    end
    
    %% Figure des données des paramètres relatif à la glace
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/'];% 
    filenamecmp=[float_list{ik} '_ice_param.pdf'];
    if exist([dir_fig,filenamecmp])
        copyfile([dir_fig,filenamecmp],[dir_tex,filenamecmp]);
    end
    
    if exist([dir_tex,filenamecmp])
        fprintf(fw1,'%s\n', '\clearpage');
        fprintf(fw1,'%s\n', ['\subsection {Ice related parameters}']);
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\begin{figure}[h!]');
        fprintf(fw1,'%s\n', '\begin{center}');
        fprintf(fw1,'%s\n', '$$');
        str=['\includegraphics[width=15cm,trim=  40 120 40 90, angle=0,clip=true]{'  filenamecmp '}'];
        fprintf(fw1,'%s\n', str);
        fprintf(fw1,'%s\n', '$$');
        str=['Float ' float_list{ik} ': Ice detection number and Transmission delayed number (upper panel), Temperature ($^{\circ}$C) in the upper 50-dbar (middle panel) and position Flag QC (lower panel), 1: good, 2: probably good, 3: probably bad, 4: bad, 8: interpolated'];
        fprintf(fw1,'%s\n', ['\caption{' str '}']);
        fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_flag}']);
        fprintf(fw1,'%s\n', '\end{center}');
        fprintf(fw1,'%s\n', '\end{figure}');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        
    end
    
    
    if length(Profflag{ik})~=0
        %%  OPTION: Figure des modifs de flags
        fprintf(fw1,'%s\n', '\clearpage');
        fprintf(fw1,'%s\n', ['\subsection {Modification of RT flags : selection of some profiles}']);
        
        %for gj=1:length(Profflag{ik})
        gj=1;
        if isstr(Profflag{ik}(gj))
            filenamecmp=[ float_list{ik} '_changeflag_' Profflag{ik}(gj,:) '.png'];
        else
            filenamecmp=[ float_list{ik} '_changeflag_' num2str(Profflag{ik}(gj)) '.png'];
        end
        
        dir_fig=[ DIR_PLOT 'change_flag/' float_list{ik} '/'];% Added By TR 28.09.20
        if exist([dir_fig,filenamecmp])
            copyfile([dir_fig,filenamecmp],[dir_tex,filenamecmp]);
        end
        
        %keyboard
        if exist([dir_tex,filenamecmp])
            
            fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            fprintf(fw1,'%s\n', '\begin{figure}[h!]');
            fprintf(fw1,'%s\n', '\begin{center}');
            fprintf(fw1,'%s\n', '$$');
            str=['\includegraphics[width=16cm,trim=  0 0 0 0, angle=0,clip=true]{'  filenamecmp '}'];
            fprintf(fw1,'%s\n', str);
            fprintf(fw1,'%s\n', '$$');
            if isstr(Profflag{ik}(gj))
                
                str=['Float ' float_list{ik} ' Comparison of Cycles ' strrep(Profflag{ik}(gj,:),'_',' and ')];
            else
                str=['Float ' float_list{ik} ' Cycle ' num2str(Profflag{ik}(gj)) '. Flag modification during DM ckeks, according to table 2'];
            end
            fprintf(fw1,'%s\n', ['\caption{' str '}']);
            fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_flag}']);
            fprintf(fw1,'%s\n', '\end{center}');
            fprintf(fw1,'%s\n', '\end{figure}');
            fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            %end
            fprintf(fw1,'%s\n', '\clearpage');
        end
        %end
    end
    
   
    if exist([dir_tex cpcor_filename])
        %% OPTION: FIGURE ANALYSE DU CPCOR
        fprintf(fw1,'%s\n', '\clearpage'); 
        fprintf(fw1,'%s\n', ['\subsection {Cpcor Analyse}']);
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\begin{figure}[h!]');
        fprintf(fw1,'%s\n', '\begin{center}');
        fprintf(fw1,'%s\n', '$$');
        
        
        fprintf(fw1,'%s\n', '$$');
        str=['\includegraphics[width=13cm,trim= 0 0 0 0, clip=true]{' cpcor_filename '}'];
        fprintf(fw1,'%s\n', str);
        % fprintf(fw1,'%s\n', '$$');
        % fprintf(fw1,'%s\n', '\end{subfigure}');
        str=['Float ' float_list{ik} '. Estimation of the optimal Cpcor and offset. Comparison with the Working Group (WG) suggested value for Cpcor.'];
        fprintf(fw1,'%s\n', ['\caption{' str '}']);
        fprintf(fw1,'%s\n', ['\label{fig4' labell '}']);
        fprintf(fw1,'%s\n', '\end{center}');
        fprintf(fw1,'%s\n', '\end{figure}');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot_dbrut%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        str=['In what follows, the salinity has been adjusted using  CPcor\_new = $' num2str(NEW_CPCOR) '$ $ dbar^{-1}$. No offset has been applied yet.']
        fprintf(fw1,'%s\n', str );
        fprintf(fw1,'%s\n', '\clearpage');
    end
    
    
    thefilename=[DIR_PLOT 'verif_profil1/' float_list{ik} '/' float_list{ik} '_T_S_bathy_ptmp0.png'];
  
    if exist(thefilename)
        %% SI DISPO: FIGURE de comparaison avec la CTD de mise a l'eau
        dir_fig=[ DIR_PLOT 'verif_profil1/' float_list{ik} '/'];% Added By TR 28.09.20
        file=[float_list{ik} '_T_S_bathy_ptmp0.png'];
        if exist([dir_fig,file])
            copyfile([dir_fig,file],[dir_tex,file]);
        end
        
        
        fprintf(fw1,'%s\n', '\clearpage');
        fprintf(fw1,'%s\n', ['\subsection {Comparison with the reference CTD cast}']);
        %fprintf(fw1,'%s\n', ['\subsubsection {Cycle 1D}']);
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\begin{figure}[h!]');
        fprintf(fw1,'%s\n', '\begin{center}');
        fprintf(fw1,'%s\n', '$$');
        
        fprintf(fw1,'%s\n', '$$');
        %str=['\includegraphics[width=16cm,trim= 50 180 60 170, clip=true]{' thefilename_pdf '}'];
        str=['\includegraphics[width=18cm,trim= 0 0 0 0, clip=true]{' file '}'];% Modified TR 28.09.2020
        fprintf(fw1,'%s\n', str);
        % fprintf(fw1,'%s\n', '$$');
        % fprintf(fw1,'%s\n', '\end{subfigure}');
        str=['Float ' float_list{ik} '. Comparison of the first descending (or ascending) argo profile  with the CTD made at float deployement. Difference is PSAL(argo) -PSAL(ref cast). '];
        fprintf(fw1,'%s\n', ['\caption{' str '}']);
        fprintf(fw1,'%s\n', ['\label{fig4' labell '}']);
        fprintf(fw1,'%s\n', '\end{center}');
        fprintf(fw1,'%s\n', '\end{figure}');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot_dbrut%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\clearpage');
        
    end
    
    %% VERIF FLAG: COMPARAISON a des profils de reference
    if length(Profref{ik})==0
        fprintf(fw1,'%s\n', '\clearpage');
        
    else
        fprintf(fw1,'%s\n', '\clearpage');
        fprintf(fw1,'%s\n', ['\subsection {Comparison to reference profiles}']);
    end
    for gj=1:length(Profref{ik})
        
        filenamecmp=''; theref='';
        if exist([DIR_PLOT 'verif_flag/' float_list{ik} '/verif_flag_' float_list{ik} '_cmpARGO_prof' num2str(Profref{ik}(gj)) '_4.png'])==2
            filenamecmp=[DIR_PLOT 'verif_flag/' float_list{ik} '/verif_flag_' float_list{ik} '_cmpARGO_prof' num2str(Profref{ik}(gj)) '_4.png'];
            theref='Argo';
            
            dir_fig=[ DIR_PLOT 'verif_flag/' float_list{ik} '/'];% Added By TR 25.09.20
            file =['verif_flag_' float_list{ik} '_cmpARGO_prof' num2str(Profref{ik}(gj)) '_4.png'];
            if exist([dir_fig,file])
                copyfile([dir_fig,file],[dir_tex,file]);
            end
            
        elseif exist([DIR_PLOT 'verif_flag/' float_list{ik} '/verif_flag_' float_list{ik} '_prof' num2str(Profref{ik}(gj)) '_4.png'])==2
            filenamecmp=[DIR_PLOT 'verif_flag/' float_list{ik} '/verif_flag_' float_list{ik} '_prof' num2str(Profref{ik}(gj)) '_4.png']
            theref='CTD';
            
            dir_fig=[ DIR_PLOT 'verif_flag/' float_list{ik} '/'];% Added By TR 25.09.20
            file =['verif_flag_' float_list{ik} '_prof' num2str(Profref{ik}(gj)) '_4.png'];
            if exist([dir_fig,file])
                copyfile([dir_fig,file],[dir_tex,file]);
            end
            
            
        end
        
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\begin{figure}[h!]');
        fprintf(fw1,'%s\n', '\begin{center}');
        fprintf(fw1,'%s\n', '$$');
        str=['\includegraphics[width=18cm,trim=  0 0 0 0, clip=true]{'  file '}'];
        fprintf(fw1,'%s\n', str);
        fprintf(fw1,'%s\n', '$$');
        str=['Float ' float_list{ik} ' Cycle ' num2str(Profref{ik}(gj)) '. The analysed Argo profile (black) is compared to the 50 nearest reference ' theref ' profiles and to two specific profiles: the nearest reference profile in time (magenta) and the nearest reference profile in space (blue). The color of reference profiles represents the year of acquisition. $\theta / S$ diagram (left panel) and a zoom on the deepest layers (right panel).'];
        fprintf(fw1,'%s\n', ['\caption{' str '}']);
        fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_refthetaS_' num2str(gj) '}']);
        fprintf(fw1,'%s\n', '\end{center}');
        fprintf(fw1,'%s\n', '\end{figure}');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        %end
        fprintf(fw1,'%s\n', '\clearpage');
    end
        %% OPTION: COMPARAISON a des profils du GDAC proches 

    if length(Profclose{ik})>0
        fprintf(fw1,'%s\n', '\clearpage');
        fprintf(fw1,'%s\n', ['\subsection {Comparison to Real Time Argo  profiles}']);
        fprintf(fw1,'%s\n', ['\begin{itemize}']);
        fprintf(fw1,'%s\n', ['\item {RT profile very close in space and time }']);
        fprintf(fw1,'%s\n', ['\end{itemize}']);
    else
        fprintf(fw1,'%s\n', '\clearpage');
    end
    for gj=1:length(Profclose{ik})
        
        if exist([DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_cycle_' num2str(Profclose{ik}(gj)) '.eps'])==2 & exist([DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_cycle_' num2str(Profclose{ik}(gj)) '.pdf'])==0
            eval(['!ps2pdf ' DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_cycle_' num2str(Profclose{ik}(gj)) '.eps ' DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_cycle_' num2str(Profclose{ik}(gj)) '.pdf'])
        end
        
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\begin{figure}[h!]');
        fprintf(fw1,'%s\n', '\begin{center}');
        fprintf(fw1,'%s\n', '$$');
        
        dir_fig=[ DIR_PLOT 'ref_database/' float_list{ik} '/'];% Added By TR 28.09.20
        file =[float_list{ik} '_cycle_' num2str(Profclose{ik}(gj)) '.png'];
        if exist([dir_fig,file])
            copyfile([dir_fig,file],[dir_tex,file]);
        end
        %str=['\includegraphics[width=16cm,trim=  50 170 50 150, clip=true]{'  DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_cycle_' num2str(Profclose{ik}(gj)) '.pdf}'];
        str=['\includegraphics[width=16cm,trim=  0 0 0 0, clip=true]{'  file '}'];
        fprintf(fw1,'%s\n', str);
        fprintf(fw1,'%s\n', '$$');
        str=['Comparison of the profile ' num2str(Profclose{ik}(gj)) ' of the float ' float_list{ik} ' with close real time Argo profiles (magenta) found within 0.1$^\circ$ in latitude, 0.2 $^\circ$ in longitude and 365 days in times. Only RT profiles with QC=1 are selected and adjusted values are used if available. Differences between the salinity of the float and the salinity of the real time Argo profiles are computed on theta level and averaged for depth below 1000db (positive values mean that the float salinities are too high)'];
        fprintf(fw1,'%s\n', ['\caption{' str '}']);
        fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_refthetaS_' num2str(gj) '}']);
        fprintf(fw1,'%s\n', '\end{center}');
        fprintf(fw1,'%s\n', '\end{figure}');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        
        
        
        %end
        fprintf(fw1,'%s\n', '\clearpage');
    end
    
    if exist([DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_and_CURR_at_theta.eps'])==2 & exist([DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_and_CURR_at_theta.pdf'])==0
        eval(['!ps2pdf ' DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_and_CURR_at_theta.eps ' DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_and_CURR_at_theta.pdf'])
    end
    
    if exist([DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_and_CURR_at_theta.png'])==2&&COMP_GDAC==1
        %fprintf(fw1,'%s\n', ['\begin{itemize}']);
        fprintf(fw1,'%s\n', ['\subsection {RT profiles in the surrounding area}']);
        %fprintf(fw1,'%s\n', ['\end{itemize}']);
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\begin{figure}[h!]');
        fprintf(fw1,'%s\n', '\begin{center}');
        fprintf(fw1,'%s\n', '$$');
        
        dir_fig=[ DIR_PLOT 'ref_database/' float_list{ik} '/' ];% Added By TR 28.09.20
        file =[float_list{ik} '_and_CURR_at_theta.png'];
        if exist([dir_fig,file])
            copyfile([dir_fig,file],[dir_tex,file]);
        end
        %str=['\includegraphics[width=16cm,trim=  50 170 50 150, clip=true]{'  DIR_PLOT 'ref_database/' float_list{ik} '/' float_list{ik} '_and_CURR_at_theta.pdf}'];
        str=['\includegraphics[width=16cm,trim=  0 0 0 0, clip=true]{'  file '}'];
        fprintf(fw1,'%s\n', str);
        fprintf(fw1,'%s\n', '$$');
        str=['Comparison of the salinity of the float ' float_list{ik} ' with  real time Argo salinities in the surrounding area at specified theta levels. Only RT data with QC=1 are selected and adjusted values are used if available. (upper panel) Map of salinities in the area. Data from float ' float_list{ik} ' are circled with magenta . (lower panel)  Salinity time series in the same area. Data from float ' float_list{ik} ' are drawn in magenta.'];
        fprintf(fw1,'%s\n', ['\caption{' str '}']);
        fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_refthetaS_' num2str(gj) '}']);
        fprintf(fw1,'%s\n', '\end{center}');
        fprintf(fw1,'%s\n', '\end{figure}');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\clearpage');
    end
    
    
    %% FIGURES OWC
    
    fprintf(fw1,'%s\n', ['\subsection {Results of the OWC method}']);
    
    fprintf(fw1,'%s\n', ['\subsubsection{Configuration}']);
    
    file =['table_config' num2str(Num_Config{ik}) '.tex'];
    if exist(['./REPORTS/' file])
        copyfile(['./REPORTS/' file],[dir_tex,file]);
    end
    fprintf(fw1,'%s\n', ['\input{table_config' num2str(Num_Config{ik}) '.tex}']);
    
    fprintf(fw1,'%s\n', ['\subsubsection{Plots}']);
    fprintf(fw1,'%s\n', ' %%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf(fw1,'%s\n','\begin{figure}[h!]');
    fprintf(fw1,'%s\n','\begin{center}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.65\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/'];% Added By TR 29.09.20
    file =[float_list{ik} '_1.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 30 150 70 100, clip=true]{' DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_1.pdf}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=clip=true]{'  file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.80\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' ];% Added By TR 29.09.20
    file =[float_list{ik} '_10.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 15 10 100 350, clip=true]{' DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_10.pdf}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{'  file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    
    str=['Float ' float_list{ik} '. Upper (left): Configuration parameters used for OWC method. Upper (rigth) : Reference  profiles used for the mapping (grey dots) are shown on the map  along with the float trajectory. Lower: Number of reference profile available within the defined spatial and temporal scales.'];
    fprintf(fw1,'%s\n', ['\caption{' str '}']);
    fprintf(fw1,'%s\n', '\label{fig1}');
    fprintf(fw1,'%s\n', '\end{center}');
    fprintf(fw1,'%s\n', '\end{figure}');
    fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf(fw1,'%s\n', ' ');
    
    fprintf(fw1,'%s\n', '\clearpage');
    
    %%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file_calib= [DIR_DATA 'float_calib/CONFIG' Num_Config{ik} '/calseries_' float_list{ik} '.mat'];
    calib.breaks='?';
    calib.max_breaks='?';
    calib.use_theta_lt='?';
    calib.use_theta_gt='?';
    calib.use_pres_lt='?';
    calib.use_pres_gt='?';
    calib.use_percent_gt='?';
    calseries=[];
    if exist(file_calib);
        calib=load(file_calib);
        calseries=calib.calseries;
        
        thefield=fieldnames(calib);
        
        for kfields=1:length(thefield)
            if isempty(calib.(thefield{kfields}))
                calib.(thefield{kfields})='[ ]';
            else
                calib.(thefield{kfields})=num2str(calib.(thefield{kfields}));
            end
        end
    end
    calseries_0='';
    calseries_1='';
    b=unique(calseries);
    f0=find(calseries==0);
    if isempty(f0)==0
        calseries_0=[num2str(f0)];
    end
    index_cal=[1:length(calseries)];
    index_cal=index_cal(~calseries==0);
    calseries=calseries(~calseries==0);
    for kcal=2:length(calseries)-1
        if calseries(kcal)~=calseries(kcal-1)
            calseries_1=[calseries_1 num2str(index_cal(kcal)) ' '];
        end
    end
    
    
    
    fprintf(fw1,'%s\n', ' %%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf(fw1,'%s\n','\begin{figure}[h!]');
    fprintf(fw1,'%s\n','\begin{center}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{1\linewidth}');
    fprintf(fw1,'%s\n', '\begin{minipage}{0.35\linewidth}');
    fprintf(fw1,'%s\n', ['\renewcommand\arraystretch{1.2}']);
    fprintf(fw1,'%s\n', '\begin{tabular}{|l|c|}');
    fprintf(fw1,'%s\n', '\hline');
    fprintf(fw1,'%s\n', [' \multicolumn{2}{|c|}{set\_calseries.m}  \\']);
    fprintf(fw1,'%s\n', '\hline');
    fprintf(fw1,'%s\n', ['breaks         & ' calib.breaks   ' \\']);
    fprintf(fw1,'%s\n', ['max\_breaks       & ' calib.max_breaks ' \\']);
    fprintf(fw1,'%s\n', [' use\_theta\_lt    & ' calib.use_theta_lt ' \\']);
    fprintf(fw1,'%s\n', [' use\_theta\_gt    & ' calib.use_theta_gt ' \\']);
    fprintf(fw1,'%s\n', [' use\_pres\_lt    & ' calib.use_pres_lt ' \\']);
    fprintf(fw1,'%s\n', [' use\_pres\_gt    & ' calib.use_pres_gt ' \\']);
    fprintf(fw1,'%s\n', ['use\_percent\_gt    & ' calib.use_percent_gt '\\']);
    
    fprintf(fw1,'%s\n', '\hline');
    fprintf(fw1,'%s\n', '\end{tabular}');
    
    fprintf(fw1,'%s\n', '\end{minipage}');
    fprintf(fw1,'%s\n', '\begin{minipage}{0.65\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' ];% Added By TR 29.09.20
    file =[float_list{ik} '_6_1_1.pdf'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 30 90 70 400, clip=true]{' DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_6_1_1.pdf}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 30 400 70 90, clip=true]{' file '}'];  % CC 17/12/2021 modified trim
    
    if exist([DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_6_1_1.pdf'])==0
        dir_fig=[ DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' ];% Added By TR 29.09.20
        file =[float_list{ik} '_6_1_2.pdf'];
        if exist([dir_fig,file])
            copyfile([dir_fig,file],[dir_tex,file]);
        end
        %str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 30 90 70 400, clip=true]{' DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_6_1_2.pdf}'];
        str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 30 400 70 90, clip=true]{'  file '}']; %CC 17/12/2021 modified trim
    end
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{minipage}');
    fprintf(fw1,'%s\n','\end{subfigure}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{1\linewidth}');
    
    fprintf(fw1,'%s\n','\begin{minipage}{0.35\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' ];% Added By TR 29.09.20
    file =[float_list{ik} '_8_1.pdf'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 350 400 50 90, clip=true]{' DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_8_1.pdf}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 350 400 50 90, clip=true]{'  file '}'];
    
    if exist([DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_8_1.pdf'])==0
        dir_fig=[ DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' ];% Added By TR 29.09.20
        file =[float_list{ik} '_8_2.pdf'];
        if exist([dir_fig,file])
            copyfile([dir_fig,file],[dir_tex,file]);
        end
        %str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 350 400 50 90, clip=true]{' DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_8_2.pdf}'];
        str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 350 400 50 90, clip=true]{'  file '}'];
    end
    fprintf(fw1,'%s\n', str);
    
    fprintf(fw1,'%s\n','\end{minipage}');
    fprintf(fw1,'%s\n','\begin{minipage}{0.65\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' ];% Added By TR 29.09.20
    file =[float_list{ik} '_9.pdf'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 30 100 70 400, clip=true]{' DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_9.pdf}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 30 100 70 400, clip=true]{'  file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{minipage}');
    
    fprintf(fw1,'%s\n','\end{subfigure}');
    str=['Float ' float_list{ik} '.  Results of the OWC method (configuration ' Num_Config{ik} '). Upper panel (right): float salinities at one $\theta$ level (blue dots) compared to mapped salinities with errors (red). Lower panel (left): The 10 $\theta$ levels (green lines) with less salinity variance along the float path that are used for computing the conductivity correction. Lower panel (right): vertically-averaged mapped salinities minus float salinities on the 10 $\theta$ levels (red) and the computed offset (green).'];
    if isempty(calseries_0)==0
        str=[str ' Omitted profiles: ' calseries_0 '.'];
    end
    if isempty(calseries_1)==0
        
        str=[str ' Time series cutted at profiles: ' calseries_1 '.'];
    end
    fprintf(fw1,'%s\n', ['\caption{' str '}']);
    fprintf(fw1,'%s\n', '\label{fig1}');
    fprintf(fw1,'%s\n', '\end{center}');
    fprintf(fw1,'%s\n', '\end{figure}');
    fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf(fw1,'%s\n', ' ');
    
    % FIGURE des comparaisons ow sur deux zones de profondeur
    %keyboard
    if exist ([DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_15.pdf'])==2 & PLOTDEPTHDEP==1
        fprintf(fw1,'%s\n', '\clearpage');
        
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\begin{figure}[h!]');
        fprintf(fw1,'%s\n', '\begin{center}');
        fprintf(fw1,'%s\n', '$$');
        
        dir_fig=[ DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' ];% Added By TR 29.09.20
        file =[float_list{ik} '_15.pdf'];
        if exist([dir_fig,file])
            copyfile([dir_fig,file],[dir_tex,file]);
        end
        %str=['\includegraphics[width=16cm,trim=  50 170 50 150, clip=true]{' DIR_PLOTOW Num_Config{ik} '/'  float_list{ik} '/' float_list{ik} '_15.pdf}'];
        str=['\includegraphics[width=16cm,trim=  50 170 50 150, clip=true]{'  file '}'];
        fprintf(fw1,'%s\n', str);
        fprintf(fw1,'%s\n', '$$');
        str=['Float ' float_list{ik} '.  Results of the OWC method (configuration ' Num_Config{ik} ') using different depths to compute calibration curve.'];
        fprintf(fw1,'%s\n', ['\caption{' str '}']);
        fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_refthetaS_' num2str(gj) '}']);
        fprintf(fw1,'%s\n', '\end{center}');
        fprintf(fw1,'%s\n', '\end{figure}');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    end
    
    %% CONCLUSION
    
    if exist(file_cor)~=0
        iil=find(findstr_tab(wmo_corr,float_list{ik}));
        fprintf(fw1,'%s\n', ['\textbf{Conclusion} ' ]);
        fprintf(fw1,'%s\n', ' ');
        fprintf(fw1,'%s\n', [ col5{iil}]);
    end
    
    fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf(fw1,'%s\n', ' ');
    fprintf(fw1,'%s\n', '\clearpage');
    fprintf(fw1,'%s\n', ['\subsection {Adjusted data}']);
    fprintf(fw1,'%s\n', ['\subsubsection {Salinity flags and correction in D files }']);
    fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf(fw1,'%s\n', '\begin{figure}[h!]');
    fprintf(fw1,'%s\n', '\begin{center}');
    fprintf(fw1,'%s\n', '$$');
    
    dir_fig=[ DIR_PLOT 'corrections/' float_list{ik} '/' ];% Added By TR 29.09.20
    file =['bilan_psal_corrections.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[width=15cm,trim=  0 40 0 0, clip=true]{'  DIR_PLOT 'corrections/' float_list{ik} '/bilan_psal_corrections.png}'];
    str=['\includegraphics[width=18cm,trim=  0 0 0 0, clip=true]{'  file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n', '$$');
    str=['Salinity correction and flags in D files (Flag 0: blue, Flag 1: green, Flag 2: yellow, Flag 3: magenta, Flag 4: red)'];
    fprintf(fw1,'%s\n', ['\caption{' str '}']);
    fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_correction}']);
    fprintf(fw1,'%s\n', '\end{center}');
    fprintf(fw1,'%s\n', '\end{figure}');
    fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    %end
    fprintf(fw1,'%s\n', '\clearpage');
    
    fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf(fw1,'%s\n', ' ');
    fprintf(fw1,'%s\n', '\clearpage');
    fprintf(fw1,'%s\n', ['\subsubsection {Sections along the float trajectory }']);
    iil=[];
    if exist(file_cor)~=0
        iil=find(findstr_tab(wmo_corr,float_list{ik}));
        fprintf(fw1,'%s\n', ['Salinity Correction applied in DM: ' col4{iil}]);
    end
    % figure des sections TPOT,PSAL,SIG0
    fprintf(fw1,'%s\n', ' %%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf(fw1,'%s\n','\begin{figure}[h!]');
    fprintf(fw1,'%s\n','\begin{center}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.48\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/' ];% Added By TR 29.09.20
    file =[float_list{ik} '_TPOT_interp_flagused_dm.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=270,origin=c,width=\textwidth,trim= 15 10 20 50, clip=true]{' DIR_PLOT 'preliminaire/' float_list{ik} '/' float_list{ik} '_TPOT_interp_flagused_dm.png}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{'  file '}'];
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.48\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/'];% Added By TR 29.09.20
    file =[float_list{ik} '_SIG0_interp_flagused_dm.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=270,origin=c,width=\textwidth,trim= 15 10 20 50, clip=true]{' DIR_PLOT 'preliminaire/' float_list{ik} '/' float_list{ik} '_PSAL_interp_flagused_dm.png}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{'  file '}'];
    
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{1\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/' ];% Added By TR 29.09.20
    file =[float_list{ik} '_PSAL_interp_flagused_dm.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=270,origin=c,width=\textwidth,trim= 15 10 20 50, clip=true]{' DIR_PLOT 'preliminaire/' float_list{ik} '/' float_list{ik} '_SIG0_interp_flagused_dm.png}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{'  file '}'];
    
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    str=['Float ' float_list{ik} '. Potential temperature, salinity and Sig0 sections along the float trajectory (adjusted data, flags used)'];
    fprintf(fw1,'%s\n', ['\caption{' str '}']);
    fprintf(fw1,'%s\n', '\label{fig1}');
    fprintf(fw1,'%s\n', '\end{center}');
    fprintf(fw1,'%s\n', '\end{figure}');
    fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf(fw1,'%s\n', ' ');
    
    fprintf(fw1,'%s\n', '\clearpage');
    fprintf(fw1,'%s\n', ['\subsubsection {Theta/S diagrams }']);
    % figure du diag OS
    fprintf(fw1,'%s\n', ' %%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf(fw1,'%s\n','\begin{figure}[h!]');
    fprintf(fw1,'%s\n','\begin{center}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.65\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/' ];% Added By TR 29.09.20
    file =[float_list{ik} '_thetaS_flagused_dm.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=270,origin=c,width=\textwidth,trim= 15 10 20 50, clip=true]{' DIR_PLOT 'preliminaire/' float_list{ik} '/' float_list{ik} '_thetaS_flagused_dm.png}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{'  file '}'];
    
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    fprintf(fw1,'%s\n','\begin{subfigure}[b]{0.65\linewidth}');
    fprintf(fw1,'%s\n','\centering');
    
    dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/' ];% Added By TR 29.09.20
    file =[float_list{ik} '_thetaS_zoom_flagused_dm.png'];
    if exist([dir_fig,file])
        copyfile([dir_fig,file],[dir_tex,file]);
    end
    %str=['\includegraphics[angle=270,origin=c,width=\textwidth,trim= 15 10 20 50, clip=true]{' DIR_PLOT 'preliminaire/' float_list{ik} '/' float_list{ik} '_thetaS_zoom_flagused_dm.png}'];
    str=['\includegraphics[angle=0,origin=c,width=\textwidth,trim= 0 0 0 0, clip=true]{'  file '}'];
    
    fprintf(fw1,'%s\n', str);
    fprintf(fw1,'%s\n','\end{subfigure}');
    str=['Float ' float_list{ik} '. Theta/S diagrams of the adjusted data, with the potential temperature referenced to 0db. Full profiles (upper panel) and zoom below 1500m (lower panel). Flags are used'];
    fprintf(fw1,'%s\n', ['\caption{' str '}']);
    fprintf(fw1,'%s\n', '\label{fig3}');
    fprintf(fw1,'%s\n', '\end{center}');
    fprintf(fw1,'%s\n', '\end{figure}');
    
    thefilename=[DIR_PLOT 'verif_profil1/' float_list{ik} '/' float_list{ik} '_T_S_zoom_ptmp0_adj.png'];
    %     thefilename_pdf=[DIR_PLOT 'verif_profil1/' float_list{ik} '/' float_list{ik} '_T_S_zoom_ptmp0_adj.pdf'];
    %
    %     if exist(thefilename_pdf,'file')==0
    %         eval(['!ps2pdf ' thefilename ' ' thefilename_pdf]);
    %     end
    
    if exist(thefilename)
        
        dir_fig=[ DIR_PLOT 'verif_profil1/' float_list{ik} '/'];% Added By TR 01.10.20
        file=[float_list{ik} '_T_S_zoom_ptmp0_adj.png'];
        if exist([dir_fig,file])
            copyfile([dir_fig,file],[dir_tex,file]);
        end
        
        
        fprintf(fw1,'%s\n', '\clearpage');
        fprintf(fw1,'%s\n', ['\subsubsection {Comparison with the reference CTD cast, adjusted profiles}']);
        %fprintf(fw1,'%s\n', ['\subsubsection {Cycle 1D}']);
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\begin{figure}[h!]');
        fprintf(fw1,'%s\n', '\begin{center}');
        fprintf(fw1,'%s\n', '$$');
        
        fprintf(fw1,'%s\n', '$$');
        str=['\includegraphics[width=18cm,trim= 0 0 0 0, clip=true]{' file '}'];
        fprintf(fw1,'%s\n', str);
        % fprintf(fw1,'%s\n', '$$');
        % fprintf(fw1,'%s\n', '\end{subfigure}');
        str=['Float ' float_list{ik} '. Comparison of the first descending (or ascending) argo profile  with the CTD made at float deployement. Difference is PSAL\_ADJUSTED(argo) -PSAL(ref cast). '];
        fprintf(fw1,'%s\n', ['\caption{' str '}']);
        fprintf(fw1,'%s\n', ['\label{fig4' labell '}']);
        fprintf(fw1,'%s\n', '\end{center}');
        fprintf(fw1,'%s\n', '\end{figure}');
        fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot_dbrut%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fw1,'%s\n', '\clearpage');
        
    end
    
    
    if length(Profrefadj{ik})==0
        fprintf(fw1,'%s\n', '\clearpage');
        
    else
        fprintf(fw1,'%s\n', '\clearpage');
        fprintf(fw1,'%s\n', ['\subsubsection {Comparison to reference profiles- Adjusted data}']);
    end
    for gj=1:length(Profrefadj{ik})
        
        filenamecmp=''; theref=''; filenamecmp1='';
        if exist([DIR_PLOT 'verif_flag/' float_list{ik} '/verif_flag_' float_list{ik} '_cmpARGO_prof' num2str(Profrefadj{ik}(gj)) '_4_adj.png'])==2
            %filenamecmp=[DIR_PLOT 'verif_flag/' float_list{ik} '/verif_flag_' float_list{ik} '_cmpARGO_prof' num2str(Profrefadj{ik}(gj)) '_4_adj.png'];
            dir_cmp=[DIR_PLOT 'verif_flag/' float_list{ik} '/'];
            fil_cmp=['verif_flag_' float_list{ik} '_cmpARGO_prof' num2str(Profrefadj{ik}(gj)) '_4_adj.png'];
            
            %filenamecmp1=[DIR_PLOT 'verif_flag/' float_list{ik} '/verif_flag_' float_list{ik} '_cmpARGO_prof' num2str(Profrefadj{ik}(gj)) '_3_adj.png'];
            dir_cmp1=[DIR_PLOT 'verif_flag/' float_list{ik} '/'];
            fil_cmp1=['verif_flag_' float_list{ik} '_cmpARGO_prof' num2str(Profrefadj{ik}(gj)) '_3_adj.png'];
            
            theref='Argo';
        elseif exist([DIR_PLOT 'verif_flag/' float_list{ik} '/verif_flag_' float_list{ik} '_prof' num2str(Profrefadj{ik}(gj)) '_4_adj.png'])==2
            %filenamecmp=[DIR_PLOT 'verif_flag/' float_list{ik} '/verif_flag_' float_list{ik} '_prof' num2str(Profrefadj{ik}(gj)) '_4_adj.png'];
            %filenamecmp1=[DIR_PLOT 'verif_flag/' float_list{ik} '/verif_flag_' float_list{ik} '_prof' num2str(Profrefadj{ik}(gj)) '_3_adj.png'];
            dir_cmp=[DIR_PLOT 'verif_flag/' float_list{ik} '/'];
            fil_cmp=['verif_flag_' float_list{ik} '_prof' num2str(Profrefadj{ik}(gj)) '_4_adj.png'];
            dir_cmp1=[DIR_PLOT 'verif_flag/' float_list{ik} '/'];
            fil_cmp1=['verif_flag_' float_list{ik} '_prof' num2str(Profrefadj{ik}(gj)) '_3_adj.png'];
            
            theref='CTD';
        end
        if exist([dir_cmp1 fil_cmp1])
            fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            fprintf(fw1,'%s\n', '\begin{figure}[h!]');
            fprintf(fw1,'%s\n', '\begin{center}');
            fprintf(fw1,'%s\n', '$$');
            
            dir_fig=dir_cmp1;% Added By TR 01.10.20
            file =fil_cmp1;
            if exist([dir_fig,file])
                copyfile([dir_fig,file],[dir_tex,file]);
            end
            
            %str=['\includegraphics[width=15cm,trim=  0 70 0 70, angle=-90,clip=true]{'  filenamecmp1 '}'];
            str=['\includegraphics[width=15cm,trim=  0 70 0 70, angle=-90,clip=true]{'  file '}'];
            fprintf(fw1,'%s\n', str);
            fprintf(fw1,'%s\n', '$$');
            str=['Float ' float_list{ik} ' Cycle ' num2str(Profrefadj{ik}(gj)) '. The analysed Argo adjusted profile (black) is compared to the 50 nearest reference ' theref ' profiles and to two specific profiles: the nearest reference profile in time (magenta) and the nearest reference profile in space (blue). The color of reference profiles represents the year of acquisition. $\theta / S$ diagram (left panel) and a zoom on the deepest layers (right panel).'];
            fprintf(fw1,'%s\n', ['\caption{' str '}']);
            fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_refthetaS_' num2str(gj) '}']);
            fprintf(fw1,'%s\n', '\end{center}');
            fprintf(fw1,'%s\n', '\end{figure}');
            fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        end
        
        if exist(filenamecmp)
            fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            fprintf(fw1,'%s\n', '\begin{figure}[h!]');
            fprintf(fw1,'%s\n', '\begin{center}');
            fprintf(fw1,'%s\n', '$$');
            
            dir_fig=dir_cmp;% Added By TR 01.10.20
            file =fil_cmp;
            if exist([dir_fig,file])
                copyfile([dir_fig,file],[dir_tex,file]);
            end
            
            %str=['\includegraphics[width=15cm,trim=  0 70 0 70, angle=-90,clip=true]{'  filenamecmp '}'];
            str=['\includegraphics[width=15cm,trim=  0 70 0 70, angle=-90,clip=true]{'  file '}'];
            fprintf(fw1,'%s\n', str);
            fprintf(fw1,'%s\n', '$$');
            str=['Float ' float_list{ik} ' Cycle ' num2str(Profrefadj{ik}(gj)) '. The analysed Argo adjusted profile (black) is compared to the 50 nearest reference ' theref ' profiles and to two specific profiles: the nearest reference profile in time (magenta) and the nearest reference profile in space (blue). The color of reference profiles represents the year of acquisition. $\theta / S$ diagram (left panel) and a zoom on the deepest layers (right panel).'];
            fprintf(fw1,'%s\n', ['\caption{' str '}']);
            fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_refthetaS_' num2str(gj) '}']);
            fprintf(fw1,'%s\n', '\end{center}');
            fprintf(fw1,'%s\n', '\end{figure}');
            fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            %end
            fprintf(fw1,'%s\n', '\clearpage');
        end
    end
    
    
end
fprintf(fw1,'%s\n', ['\end{document}']);

fclose(fw1);




%% genere le fichier config pour les corrections
theref_file= ['../table_versionOW_versionbase.txt'];
fr=fopen(theref_file);
versionbase=textscan(fr,'%s\n','Delimiter',',','CommentStyle','#');
%[num_ligne,dirow,verow,verctd,verargo]=get_txtfile_col(theref_file,',');
[num_ligne,verow,verctd,verargo]=get_txtfile_col(theref_file,',');
%iver=findstr_tab(dirow,CONF.VERSION_OWC);
iver=1;
fclose(fr);
if isempty(iver)
    error(['Your OWC version is not listed in this file: ' theref_file] )
else
    for ik=1:length(float_list)
        theconfig_file = [DIR_CODES '/CORRECTIONS/paramlog/config_' float_list{ik} '.txt'];
		if exist([DIR_CODES '/CORRECTIONS/paramlog/'])==0
			mkdir([DIR_CODES '/CORRECTIONS/paramlog/'])
		end
        theconfig_owc  = [DIR_CODES '/LPO_CODES_ATLN_NEW/ow_config/ow_config_' Num_Config{ik} '.txt'];
        ow_conf=load_configuration(theconfig_owc);
        
        switch ow_conf.CONFIG_WMO_BOXES
            case {'wmo_boxes_ctd.mat'}
                refdatabase=verctd{iver};
            case {'wmo_boxes_argo.mat'}
                refdatabase=verargo{iver};
            case {'wmo_boxes_ctdandargo.mat'}
                refdatabase=[verctd{iver} ' & ' verargo{iver}]
        end
        
        %theref_file= ['./REPORTS/ref_config' Num_Config{ik} '.txt'];
        fw1=fopen(theconfig_file,'w');
        % fr=fopen(theref_file);
        % refdatabase=textscan(fr,'%s','Delimiter',',','CommentStyle','#');
        
        fprintf(fw1,'%s\n', '% MAIN_write_dmqc_files configuration file');
        fprintf(fw1,'%s\n', '%%% IN/OUT DIRECTORY');
        fprintf(fw1,'%s\n', '   ');
        
        fprintf(fw1,'%s\n', '% INPUT NETCDF files DIRECTORY ( MAIN_write_dmqc_files will look in DIR_FTP/$flt_name$/profiles/)');
        fprintf(fw1,'%s\n', ['DIR_FTP=' DIR_FTP tabdac{ik} '/']);
        fprintf(fw1,'%s\n', '   ');
        
        fprintf(fw1,'%s\n', '% CALIBRATION files from OW: cal_$float_name$.mat ( MAIN_write_dmqc_files will look in DIR_OW/float_calib/)');
        fprintf(fw1,'%s\n', ['DIR_DATA=' DIR_DATA 'float_calib/CONFIG' Num_Config{ik} '/']);
        fprintf(fw1,'%s\n', '   ');
        
        fprintf(fw1,'%s\n', '% output files with DMQC corrections are put in this directory (MAIN_write_dmqc_files will put files in DIR_OUT/$flt_name$/profiles/)');
        fprintf(fw1,'%s\n', ['DIR_DM_FILES=' DIR_DMQC tabdac{ik} '/' ])
        fprintf(fw1,'%s\n', '   ');
        
        fprintf(fw1,'%s\n', '% plot from final checks');
        fprintf(fw1,'%s\n', ['DIR_PLOT=' DIR_PLOT '/corrections/' float_list{ik} '/' ])
        fprintf(fw1,'%s\n', '   ');
        
        
        fprintf(fw1,'%s\n', '%%% INFORMATIONS ON OW METHOD  (default values used in calibration comments)');
        fprintf(fw1,'%s\n', '   ');
        fprintf(fw1,'%s\n', ['VERSION= ' verow{iver}]);
        %fprintf(fw1,'%s\n', ['BASEREF= ' refdatabase{1}{1}]);
        fprintf(fw1,'%s\n', ['BASEREF= ' refdatabase]);
        fprintf(fw1,'%s\n', 'REPORT=');
        %%% INFORMATIONS ON  DMQC OPERATOR (to be written in global attributes :comment_dmqc_operator = "PRIMARY | OPERATOR_ORCID_ID | OPERATOR_NAME, OPERATOR_INSTITUTION") ;
        %----   -----------------------------
        fprintf(fw1,'%s\n', ['OPERATOR_ORCID_ID= ' CONF.OPERATOR_ORCID_ID]);
        fprintf(fw1,'%s\n', ['OPERATOR_NAME= ' CONF.OPERATOR_NAME]);
        fprintf(fw1,'%s\n', ['OPERATOR_INSTITUTION= ' CONF.OPERATOR_INSTITUTION]);
        fclose(fw1);
        
        generate_conf_table(str2num(Num_Config{ik}),refdatabase)
        % 		if ~exist([DIR_CODES '/DOC/OVERLEAF/table_config_' Num_Config{ik} '.pdf'])
        % 		   disp ('You should run this code:')
        % 		   disp(['latex ' DIR_CODES '/DOC/OVERLEAF/table_config' Num_Config{ik} '.tex'])
        % 		 end
        
    end
end

disp ('You should run this code:')
disp(['latex ' DIR_CODES '/DOC/OVERLEAF/ebauche_rapport.tex'])

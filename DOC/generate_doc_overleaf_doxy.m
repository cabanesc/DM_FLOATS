% -========================================================
%   USAGE : generate_doc_overleaf_doxy(tabfloat,tabdac)
%   PURPOSE : generate a tex file to repport DMQC of DOXY parameter
% -----------------------------------
%   INPUT :
%    tabfloat  (char or cell of chars -size n_floatsx1)    e.g. '6900258' or {'6900258', '3901954'}
%    tabdac    (char or cell of chars -size n_floatsx1)    e.g. 'coriolis' or {'coriolis', 'bodc'}
%
%   OPTIONNAL INPUT :
%
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created (2016) ccabanes
%
%   CALLED SUBROUTINES:
% -------------------------------------
%
% ========================================================% 
function generate_doc_overleaf_doxy(tabfloat,tabdac)


if iscell(tabfloat)==0;tabfloat=cellstr(tabfloat);end
if iscell(tabdac)==0;tabdac=cellstr(tabdac);end
if length(tabfloat)>1&length(tabdac)==1
    tabdac=repmat(tabdac,1,length(tabfloat));
end
float_list=tabfloat;

%init_path
C = load_configuration('config.txt');
DIR_FTP=C.DIR_FTP;
DIR_PLOT=C.DIR_PLOT;
DIR_DATA=C.DIR_DATA;
DIR_PLOTOW=[C.DIR_DATA 'float_plots/CONFIG'];
DIR_DMQC=C.DIR_DM_FILES;
DIR_CODES=C.DIR_CODES;

% OVERLEAF est un repertoire temporaire où on stoke les figures et le fichier tex du rapport
% cela facilite la compilation si on le fait via  le site overleaf
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

% Fichier Tex du rapport
fw1=fopen(['./OVERLEAF/doxy_report.tex'],'w');

% Entete fichier tex /package /format
%============================================================
fprintf(fw1,'%s\n', ['\documentclass[11pt,titlepage]{article}']);
fprintf(fw1,'%s\n', ['\usepackage[latin1]{inputenc}']);
fprintf(fw1,'%s\n', ['\usepackage[dvips,final]{graphicx}']);
fprintf(fw1,'%s\n', ['\usepackage{color}']);
fprintf(fw1,'%s\n', ['\usepackage{tabularx}']);
fprintf(fw1,'%s\n', ['\usepackage{longtable}'])
fprintf(fw1,'%s\n', ['\usepackage{float}']);
fprintf(fw1,'%s\n', ['\usepackage{calc}']);
fprintf(fw1,'%s\n', ['\usepackage{soul}']);
fprintf(fw1,'%s\n', ['\usepackage[small,hang]{caption}']);
fprintf(fw1,'%s\n', ['\usepackage{subcaption}']);
fprintf(fw1,'%s\n', ['\usepackage[colorlinks={true}]{hyperref}']);
fprintf(fw1,'%s\n', [' ']);
fprintf(fw1,'%s\n', ['\makeatletter']);
fprintf(fw1,'%s\n', ['\def\nobreakhline{\multispan\LT@cols\unskip\leaders\hrule\@height\arrayrulewidth\hfill\cr\noalign{\penalty10000}}']);
fprintf(fw1,'%s\n', ['\makeatother']);
fprintf(fw1,'%s\n', [' ']);
fprintf(fw1,'%s\n', ['\setlength{\textwidth}{430pt}']);
fprintf(fw1,'%s\n', ['\setlength{\oddsidemargin}{0pt}']);
fprintf(fw1,'%s\n', ['\setlength{\evensidemargin}{0pt}']);
fprintf(fw1,'%s\n', ['\setlength{\marginparwidth}{0pt}']);
fprintf(fw1,'%s\n', ['\setlength{\headsep}{35pt}']);
fprintf(fw1,'%s\n', ['\setlength{\topmargin}{-15pt}']);
fprintf(fw1,'%s\n', ['\setlength{\textheight}{650pt}']);
fprintf(fw1,'%s\n', ['\newcommand{\dg}{$^\circ$}']);
fprintf(fw1,'%s\n', ['\newcommand{\Sv}{$\times 10^6$ m$^3$.s$^{-1}$}']);
fprintf(fw1,'%s\n', ['\newcommand{\pw}{$\times 10^{15}$ W}']);
fprintf(fw1,'%s\n', ['\newcommand{\cmps}{cm.s$^{-1}$}']);
fprintf(fw1,'%s\n', ['\newcommand{\mps}{m.s$^{-1}$}']);
fprintf(fw1,'%s\n', ['\newcommand{\mcps}{m$^2$.s$^{-1}$}']);
fprintf(fw1,'%s\n', ['\newcommand{\upv}{$\times 10^{-11}$ m$^{-1}$.s$^{-1}$}']);
fprintf(fw1,'%s\n', ['\newcommand{\usig}{kg.m$^{-3}$}']);


fprintf(fw1,'%s\n', ['\begin{document}']);

fprintf(fw1,'%s\n', ['\pagestyle{plain}']);
fprintf(fw1,'%s\n', ['\pagenumbering{arabic}']);

% PAGE de TITRE 
%============================================================
TITRE=' CECI est un titre';
SOUS_TITRE='CECI est un sous-titre';
AUTHORS='xxx and yxy'
fprintf(fw1,'%s\n', ['\begin{titlepage}']);


fprintf(fw1,'%s\n', ['\begin{center}']);
fprintf(fw1,'%s\n', ['{']);
fprintf(fw1,'%s\n', ['\newfont{\gtitre}{cmssbx10 at 16pt}']);
fprintf(fw1,'%s\n', ['\newfont{\mtitre}{cmss10 at 14pt}']);

fprintf(fw1,'%s\n', ['\begin{tabularx}{\textwidth}{lXr}']);

fprintf(fw1,'%s\n', ['\vspace*{2cm} & & \\']);
fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['\vspace*{0.5cm} & & \\']);
fprintf(fw1,'%s\n', ['\multicolumn{3}{c}{\gtitre ' TITRE  ' } \\']);
fprintf(fw1,'%s\n', ['\vspace*{-0.2cm} & & \\']);
fprintf(fw1,'%s\n', ['\multicolumn{3}{c}{\gtitre  ' SOUS_TITRE '} \\']);
fprintf(fw1,'%s\n', ['\vspace*{0.5cm} & & \\']);
fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['\vspace*{0.5cm} & & \\']);
fprintf(fw1,'%s\n', ['\multicolumn{3}{c}{ ' AUTHORS ' }\\ ']);
fprintf(fw1,'%s\n', ['\vspace*{0.5cm} & & \\']);
fprintf(fw1,'%s\n', ['\multicolumn{3}{c}{ SO ARGO - LOPS report - Update \today}\\  ']);
fprintf(fw1,'%s\n', ['\end{tabularx}']);
fprintf(fw1,'%s\n', [' ']);
fprintf(fw1,'%s\n', ['\vspace*{2cm} ']);
fprintf(fw1,'%s\n', [' \textbf{Summary}']);
fprintf(fw1,'%s\n', ['}']);
fprintf(fw1,'%s\n', ['\end{center}']);

fprintf(fw1,'%s\n', ['\end{titlepage}']);

%  TABLE MATIERE
%============================================================
fprintf(fw1,'%s\n', ['\tableofcontents']);
fprintf(fw1,'%s\n', ['\clearpage']);
fprintf(fw1,'%s\n', ['\setcounter{section}{0}']);

% EXEMPLE:  section 
%===================
fprintf(fw1,'%s\n', ['\section{Exemple section}']);

% EXEMPLE : SUBSECTION
fprintf(fw1,'%s\n', ['\subsection{Exemple sous-section}']);
%


%----------------------------------------------------------------------------
%  ecriture de la TABLE  : Correction des flags
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
NcVar.doxy_qc.name='DOXY_QC';
NcVar.pres.name='PRES';

disp('Table 2: recherche des informations dans l''history')
% Recherche des informations dans l'history
for ik=1:length(float_list)
    %[Co,Dim]=create_multi_from_mono(DIR_FTP,float_list{ik},tabdac{ik},'CR','Primary sampling',NcVar);
    [file_list]=select_float_files_on_ftp(float_list{ik},tabdac{ik},DIR_FTP,'B');
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
% fin de la table sur la correction des flags
%----------------------------------------------------------------------------

% %% EXEMPLE : ajout d'une figure
% %===============================

% %%  filenamecmp: emplacement de la figure + nom de la figure
% %%  --------------------------------------------------------
% dir_fig=[ DIR_PLOT 'preliminaire/' float_list{ik} '/'];
% filenamecmp=[float_list{ik} '_surface_pres.pdf'];

%  %%  on copie la figure dans le repertoire temporaire ./OVERLEAF
% %%  --------------------------------------------------------

% if exist([dir_fig,filenamecmp])
	% copyfile([dir_fig,filenamecmp],[dir_tex,filenamecmp]);
% end

% %%% on ecrit dans le fichier tex 
% %%  --------------------------------------------------------
   
%fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); 
%fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%fprintf(fw1,'%s\n', '\begin{figure}[h!]');
%fprintf(fw1,'%s\n', '\begin{center}');
%fprintf(fw1,'%s\n', '$$');
%str=['\includegraphics[width=15cm,trim=  20 80 20 80, angle=0,clip=true]{'  filenamecmp '}'];
%fprintf(fw1,'%s\n', str);
%fprintf(fw1,'%s\n', '$$');
%str=['Float ' float_list{ik} ': Some technical data as read in the technical file'];  % Legende de la figure
%fprintf(fw1,'%s\n', ['\caption{' str '}']);
%fprintf(fw1,'%s\n', ['\label{fig' float_list{ik} '_flag}']);
%fprintf(fw1,'%s\n', '\end{center}');
%fprintf(fw1,'%s\n', '\end{figure}');
%fprintf(fw1,'%s\n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');




fprintf(fw1,'%s\n', ['\end{document}']);

fclose(fw1);

disp ('You should run this code:')
disp(['latex ' DIR_CODES '/DOC/OVERLEAF/doxy_report.tex'])

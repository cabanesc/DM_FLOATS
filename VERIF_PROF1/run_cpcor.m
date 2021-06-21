 %==========================================================================
% Auteur: DAVID Nicolas - IFREMER/ENSIETA
% Date: 09/08/2007
% Objectif: Comparer le profil 0A du flotteur ARGO aux donnees de la
% campagne oceanographique
% Sous-routine principale: routine_profil1.m, graphique_profil1.m,
% difference_profil1.m
% Modification: 
%==========================================================================
%close all

addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/Lib_forCathy/ObsInSitu/Lib_Argo/')
addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/Lib_forCathy/ObsInSitu/Lib_Argo/RWnetcdf/R2008b')
addpath(genpath('/home1/triagoz/matlab/outils_matlab/seawater/gsw_matlab_v3_04_TR'))

lo_system_configuration = load_configuration( 'config_bocats.txt' );
%lo_system_configuration = load_configuration( 'config_ovide18.txt' );
 %lo_system_configuration = load_configuration( 'config_rrex17.txt' );
% lo_system_configuration = load_configuration( 'config_bocats.txt' );
 %lo_system_configuration = load_configuration( 'config_rrex15.txt' );
%lo_system_configuration = load_configuration( 'config_geovide.txt' );
%lo_system_configuration = load_configuration( 'config_jbsalle18.txt' )
%lo_system_configuration = load_configuration( 'config_pirata18.txt' );

%--------------------------------------------------------------------------
% Entree des parametres
%--------------------------------------------------------------------------

if exist('auto','var') == 0 % Cas d'execution automatique via routine
%      camp_name = input('Quel est le nom de la campagne oceanographique ? ','s');
%      campagne = input('Quel est le nom du fichier .mat de la campagne (ovid04/ovid06/...) ?','s');  
%      float = input('Quel(s) flotteur(s) ? ');        
%      palier = input('Quel est le palier ? ');
%      pas_surface = input('Quel est le pas en surface ? ');
%      pas_fond = input('Quel est le pas en-dessous du palier ? ');
    
    %float = [6901763 6901601 6902818 6902881 6902882] % ovide18
    %float = [6901603 6902810 6902811 6902812 6902819] % rrex17
    %float = [6901603 6902810   6902819] % rrex17
    %float = [6901760 6901762];  % bocats
    %float = [ 6902814]; 
    %float = [6901602 6901758 6901759 6901757] ; %rrex15
    %float = [ 6901758 6901759 6901757] ; %rrex15 cycle 1A
    %float = [6901631 6901632] ; % geovide
    %float =[3902131 3902132]; % pirata18
    %float =[3902129 6902813 6902814]; % jbsalle18
    % float =[6901246 6901248 3902126 3902127] % pedro velez floats
    %mindepth ={'2500','2500','2000','1500'} % rrex profil 1A
    float = {'3902127'};
    mindepth=repmat({'2000'},length(float),1);
    palier = 1500;
    pas_surface = 10;
    pas_fond = 20;
    %nocycl = input('Numero de profil  ','s');
    %float = [6902814];
    %nocycl={'5','6'}; % geovide
    nocycl=repmat({'2'},length(float),1);
    
end

[num_ligne,float,nocycl,mindepth,config_file]=get_txtfile_col([lo_system_configuration.PROF1_DIRECTORY 'liste_deep_floats_for_cpcor.csv'],',');

% 6901631, données <100m
% 6901632, données <0m
%--------------------------------------------------------------------------
% Traitement
%--------------------------------------------------------------------------
fw=fopen([lo_system_configuration.PROF1_DIRECTORY '/temp_table.csv'],'w')


%for i = 1:size(float,2)
for i = 34:34
    lo_system_configuration = load_configuration( strtrim(config_file{i}));

    %for i = 1:1
    %flt_name = num2str(float(i));
    flt_name = float{i};
    disp([datestr(now) ' Working on float ' flt_name])
    lo_system_configuration.ZOOM=mindepth{i};
	lo_system_configuration.ZOOM='2000';
    %--------------------------------------------------------------------------
    % Préparation du répertoire de sauvegarde
    %--------------------------------------------------------------------------
    dir_enregistre = [lo_system_configuration.PROF1_DIRECTORY  flt_name];
    
    if ~exist(dir_enregistre,'dir')
        mkdir(dir_enregistre);
    end

    routine_profil1_cpcor_prescor(fw,nocycl{i},flt_name,lo_system_configuration,palier,pas_surface,pas_fond);
    
end

%clear all
%close all

%=======================================================================
% mlt_lect_tab      - Lecture de la table de parametres presente dans data/wk_service/
%                     Appelee dans Coriolis.
%
%
%  EXAMPLE:
%  function [imes, NIV_STD, dat_est, ampl_AO, nb_niv, INIV_CAL, PARAM_ANO, ...
%  LS_MS_Sbgr_weight, fic_ATLAS_conv, fic_ATLAS_AO, FIC_AOB_prefx, str_doc_ini] ...
%     = mlt_lect_tab(itr);
%
%  INPUT PARAMETERS :
%      itr = 1 : lecture des paramètres concernant 
%              uniquement la conversion des Multiprofils
%      itr = 2 : lecture des paramètres concernant 
%              la conversion des Multiprofils et 
%              l'analyse objective
%
%  OUTPUT PARAMETERS :
%       imes    : mis à 1 si problème de lecture
%               de la table des paramètres - Arret du traitement
%       NIV_STD : niveaux grille (niveaux d interpolation)
%       dat_est : date d estimation (AO)
%       ampl_AO : nb de jours pris en compte dans l AO (avant et apres)
%       nb_niv : nd de niveaux d AO
%       INIV_CAL : indices tq NIV_STD == Niveaux AO (NIV_CAL)
%       PARAM_ANO : parametre a analyser, ex :'TEMP'
%       LS_MS_Sbgr_weight : Poids des amplitudes covariances Longscale,mesoscale et sous_grille
%       fic_ATLAS_conv : atlas de reference pour la conversion mlt (fournit la densite)
%       fic_ATLAS_AO : climato de reference pour les AO
%       FIC_AOB_prefx : prefixe du nom du fichier .nc resultat de l AO
%       str_doc_ini : documentation a ecrire dans les fichiers .nc resultats
%
%
%  AUTHORS : C. Lagadec, E. Autret
%
%  HISTORIC DEVELOPMENT
%      11/00 - C. Lagadec - creation Version V1.00
%      03/02 - E.Autret - Version V2.00
%
%  NOTES:
%    format de la table
%    ------------------
%       des lignes 1 à 5 : zones multiprofils
%   ligne 1 : nom du fichier ATLAS pour conversion mlt profile
%   ligne 2 : répertoire référence (contenant les fichiers 
%              d'isocontours, de bathymétrie et Atlas)
%   ligne 3 : répertoire d'importation
%   ligne 4 : répertoire de création des multiprofils
%   ligne 5 : limites de zones à récupérer (zone_RECUP)
%   ligne 6 : paramètre de référence (DEPH)
%   ligne 7 : niveaux d interpolation (NIV_STD)
%
%       des lignes 8 à 20 : Analyse Objective 
%   ligne 8 : climato de reference pour l AO (analyses sur anomalies /t climato)
%   ligne 9 : répertoire de resultats d AO
%   ligne 10 : pas en latitude-longitude pour AO
%            si il n est pas specifie, on prend le pas de fic Atlas AO
%
%========================================================================

function[imes, NIV_STD, dat_est, ampl_AO, nb_niv, INIV_CAL, PARAM_ANO, LS_MS_Sbgr_weight, fic_ATLAS_conv, fic_ATLAS_AO, FIC_AOB_prefx, str_doc_ini] = mlt_lect_tab(itr,fic_tab_AO,PARAM)


mlt_globalrep;

 NIV_STD     = [];
 dat_est     = [];
 ampl_AO     = [];
 nb_niv      = [];
 INIV_CAL    = [];
 PARAM_ANO   = [];
 LS_MS_Sbgr_weight   = [];
 fic_ATLAS_conv    = [];
 fic_ATLAS_AO      = [];
 FIC_AOB_prefx     = [];
 str_doc_ini     = [];


tab_label = strvcat('fic Atlas           ', ...
                    'rep_reference       ', ...
                    'rep_importation     ', ...
                    'rep_multiprofils    ', ...
                    'limites recup       ', ...
                    'param ref           ', ...
                    'niveaux grille      ', ...
                    'fic Atlas AO        ', ...
                    'rep_anaobj          ', ...
                    'pas                 ', ...
                    'niveaux a analyser  ', ...
                    'niveaux a tracer    ', ...  
                    'Poids LS, MS, Sbgr  ', ...  
                    'dat_est             ', ...
                    'ampl_AO             ', ...
                    'title               ', ...
                    'Experiment name     ', ...
                    'Project name        ', ...  
                    'Data manager        ', ...  
                    'Geographic Area     ', ... 
                    'fin                 ');

mess1 = 'Table mal constituée !';
mess_arret = 'Arrêt du traitement';


rep_mlt_AO  = [];
 

itab = 0;

imess_choixtab = 0;		% mis à 1 si pas de table sélectionnée
                                % dans AO_choix_tab

imes = 0;		        % mis à 1 pour ne pas continuer 
                              % le traitement dans coriolis
			      

			      
if (nargin == 1)
  [fic_tab_AO, imess_choixtab] = AO_choix_tab;
  PARAM=[];
else
  fic_tab_AO=[pwd '/wk_service/' fic_tab_AO];
end

fprintf(fic_tab_AO);


if  imess_choixtab == 0

          ftab_AO = fopen(fic_tab_AO,'r');



% ligne 1 : fichier Atlas
% -----------------------

         tab_AO = fgetl(ftab_AO);
         label_AO = tab_AO(1:20);
         if   label_AO ~= tab_label(1,:)
                  h= warndlg(mess1,mess_arret);
                  waitfor(h);
                  break;
         end
	ltab = length(tab_AO);
	fic_ATLAS_conv = tab_AO(21:ltab);

	if ~exist(fic_ATLAS_conv,'file'),
         	h= warndlg([fic_ATLAS_conv ' : Fichier inexistant'],mess_arret);
         	waitfor(h);
        	 break;
	end



% ligne 2 : répertoire référence
% ------------------------------

	tab_AO = fgetl(ftab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(2,:)
         	h= warndlg(mess1,mess_arret);
         	waitfor(h);
        	 break;
	end
	ltab = length(tab_AO);

	rep_iso = tab_AO(21:ltab);
	l = length(rep_iso);
	if  rep_iso(l:l) ~= '/'
         	rep_iso = [rep_iso '/'];
	end


	if  ~exist(rep_iso,'dir')  
    		h = warndlg([rep_iso ,' : le répertoire des fichiers ISOCONTOURS 	n''existe pas.'],'Arrêt du traitement');
   		 break;
	end
	
% ligne 3 : repertoire d'origine pour les fichiers a importer
%------------------------------------------------------------

	tab_AO = fgetl(ftab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(3,:)
         	h= warndlg(mess1,mess_arret);
         	waitfor(h);
         	break;
	end
	ltab = length(tab_AO);

        rep_camp     = deblank(tab_AO(21:ltab));
        lr = length(rep_camp);
        if rep_camp(lr) ~= '/'
           rep_camp = [rep_camp '/'];
        end;

% ligne 4 :  repertoire  d'écriture des fichiers interpolés
%----------------------------------------------------------

        rep_ANE     = [pwd '/'];

	tab_AO = fgetl(ftab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(4,:)
         	h= warndlg(mess1,mess_arret);
         	waitfor(h);
         	break;
	end

	ltab = length(tab_AO);

	rep_mlt_suf  = deblank(tab_AO(21:ltab));

	if rep_mlt_suf(length(rep_mlt_suf)) ~= '/'
        	rep_mlt_suf = [rep_mlt_suf '/'];
	end

	rep_mlt = [rep_ANE rep_mlt_suf];
	if ~exist(rep_mlt,'dir'),
             h = warndlg('Le répertoire d''écriture des fichiers Multistations n''existe pas. Création automatique !','Attention !');
             mkdir(rep_ANE,rep_mlt_suf);
	end

% création (ou non) des répertoires implicites wk_...

	mlt_initrep;



% ligne 5 : limites de zone à récupérer
% -------------------------------------

	tab_AO = fgetl(ftab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(5,:)
         	h= warndlg(mess1,mess_arret);
         	waitfor(h);
         	break;
	end
	ltab = length(tab_AO);

	zone_RECUP   = tab_AO(21:ltab);
	zone_RECUP   = str2num(zone_RECUP)';


% ligne 6 : parametre de reference (param ref)
% --------------------------------------------

tab_AO = fgetl(ftab_AO);
label_AO = tab_AO(1:20);
if   label_AO ~= tab_label(6,:)
         h= warndlg(mess1,mess_arret);
         waitfor(h);
         break;
end

% ne sert que pour l'écriture du paramètre
% de référence du fichier Résultat de l'analyse objective

MLT_PAREF= tab_AO(21:24);


% ligne 7 : niveaux de la grille multiprofils
% -------------------------------------------

	tab_AO   = fgetl(ftab_AO);
	ltab = length(tab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(7,:)
        	h= warndlg(mess1,mess_arret);
        	 waitfor(h);
        	break;
	end 

	NIV_STD = str2num(tab_AO(21:ltab))';

	nb_std = length(NIV_STD);

	INIV_STD = [1:nb_std];






% ===========================================================
%
% lecture des paramètres Analyse Objective
%
% ===========================================================


if  itr == 2 

% ligne 8 : fichier Atlas AO
% -----------------------

         tab_AO = fgetl(ftab_AO);
         label_AO = tab_AO(1:20);
         if   label_AO ~= tab_label(1,:)
                  h= warndlg(mess1,mess_arret);
                  waitfor(h);
                  break;
         end
	ltab = length(tab_AO);
	fic_ATLAS_AO = tab_AO(21:ltab);

	if ~exist(fic_ATLAS_AO,'file'),
         	h= warndlg([fic_ATLAS_AO ' : Fichier inexistant'],mess_arret);
         	waitfor(h);
        	 break;
	end
	str_doc_ini= strvcat(str_doc_ini, tab_AO(21:ltab));

% ligne 9 :  repertoire  d'écriture des fichiers AO
%----------------------------------------------------------

    rep_ANE     = [pwd '/'];

	tab_AO = fgetl(ftab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(4,:)
         	h= warndlg(mess1,mess_arret);
         	waitfor(h);
         	break;
	end

	ltab = length(tab_AO);

	rep_ao_suf  = deblank(tab_AO(21:ltab));

	if rep_ao_suf(length(rep_ao_suf)) ~= '/'
        	rep_ao_suf = [rep_ao_suf '/'];
	end

	rep_anaobj = [rep_ANE rep_ao_suf];
	if ~exist(rep_anaobj,'dir'),
             h = warndlg('Le répertoire d''écriture des fichiers AO n''existe pas. Création automatique !','Attention !');
             mkdir(rep_ANE,rep_ao_suf);
	end

    


% ligne 10 :  pas (latitude et longitude) 
% --------------------------------------

	tab_AO = fgetl(ftab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(8,:)
         	h= warndlg(mess1,mess_arret);
         	waitfor(h);
         	break;
	end
    
	lat_pas  = str2num(tab_AO(21:23));
	lon_pas  = str2num(tab_AO(31:33));
    

% ligne 11 : niveaux de calcul AO 
% ------------------------------

	tab_AO   = fgetl(ftab_AO);
	ltab = length(tab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(9,:)
         	h= warndlg(mess1,mess_arret);
         	waitfor(h);
         	break;
	end 

	NIV_CAL = str2num(tab_AO(21:ltab))';

	nb_niv = length(NIV_CAL);

	for kk=1:nb_niv
  		INIV_CAL  (kk)= find(NIV_STD==NIV_CAL(kk));
	end
% 
% % ligne 12 : niveaux a tracer 
% % ---------------------------
% 
% 	tab_AO   = fgetl(ftab_AO);
% 	ltab = length(tab_AO);
% 	label_AO = tab_AO(1:20);
% 	if   label_AO ~= tab_label(10,:)
%         	h= warndlg(mess1,mess_arret);
%         	waitfor(h);
%         	break;
% 	end 
% 
% 	NIV_PRN = str2num(tab_AO(21:ltab));
% 	for kk=1:length(NIV_PRN);
%   		INIV_PRN (kk) = find(NIV_CAL==NIV_PRN(kk));
% 	end
% 

% ligne 12 : paramètre à tracer (soit TEMP,PSAL,SIGI) 
% ------------------------------------------------------

	tab_AO    = fgetl(ftab_AO);
	label_AO  = tab_AO(1:20);
	ltab = length(tab_AO);
	if   label_AO ~= tab_label(11,:)
        	 h= warndlg(mess1,mess_arret);
         	waitfor(h);
        	 break;
	end
	if isempty(PARAM)
	  PARAM_ANO = deblank(tab_AO(21:ltab));
	else
	  PARAM_ANO = PARAM;
	end
	
	
    
% ligne 13 : Poids echelles 
% ---------------------------

	tab_AO   = fgetl(ftab_AO);
	ltab = length(tab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(10,:)
        	h= warndlg(mess1,mess_arret);
        	waitfor(h);
        	break;
	end 

	LS_MS_Sbgr_weight = str2num(tab_AO(21:ltab));



% ligne 14 : date d estimation (dat_est)
% ------------------------------------

	tab_AO   = fgetl(ftab_AO);
	ltab = length(tab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(12,:)
            h= warndlg(mess1,mess_arret);
            waitfor(h);
            break;
	end

if ltab>20
    dat_est=[str2num(tab_AO(21:22));str2num(tab_AO(24:25));str2num(tab_AO(27:30))];
end
  


% ligne 15 : nombre d estimations (nb_est)
% ------------------------------------
%	tab_AO   = fgetl(ftab_AO);
%	ltab = length(tab_AO);
%	label_AO = tab_AO(1:20);
%	if   label_AO ~= tab_label(12,:)
%            h= warndlg(mess1,mess_arret);
%            waitfor(h);
%            break;
%	end

%    nb_est= str2num(tab_AO(21:ltab));



% ligne 15 : amplitude en nb jours pour Analyse Objective (ampl_AO)
% -----------------------------------------------------------------

	tab_AO   = fgetl(ftab_AO);
	ltab = length(tab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(13,:)
            h= warndlg(mess1,mess_arret);
            waitfor(h);
            break;
	end

	ampl_AO= str2num(tab_AO(21:ltab));


 	


%  ligne 16 : title
% -----------------

	tab_AO   = fgetl(ftab_AO);
	ltab = length(tab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(15,:)
            h= warndlg(mess1,mess_arret);
            waitfor(h);
            break;
	end

	str_doc_ini= strvcat(str_doc_ini, tab_AO(21:ltab));
	
	str_doc_ini = strvcat(str_doc_ini,'Software [AO-V2.0]');
    
    str_doc_ini = strvcat(str_doc_ini,'Estimate',' ');
    
    
    
    %prefixe nom du fichier AO
    %-------------------------
    ih = findstr(tab_AO(21:ltab),'weekly');
    im = findstr(tab_AO(21:ltab),'monthly');
    if ~isempty(ih)
        FIC_AOB_prefx = [rep_anaobj PARAM_ANO '/provh']
    end
    if ~isempty(im)
        FIC_AOB_prefx = [rep_anaobj PARAM_ANO '/provm']
    end
    if ~isempty(ih)
        str_doc_ini = strvcat(str_doc_ini,'weekly');
    end
    if ~isempty(im)
        str_doc_ini = strvcat(str_doc_ini,'monthly');
    end
    
    
 %  ligne 17 : Experiment name
% ------------------------
	tab_AO   = fgetl(ftab_AO);
	ltab = length(tab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(18,:)
            h= warndlg(mess1,mess_arret);
            waitfor(h);
            break;
	end
	
	str_doc_ini= strvcat(str_doc_ini, tab_AO(21:ltab));

    
    
    
    
%  ligne 18 : Project name
% ------------------------
	tab_AO   = fgetl(ftab_AO);
	ltab = length(tab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(18,:)
            h= warndlg(mess1,mess_arret);
            waitfor(h);
            break;
	end
	
	str_doc_ini= strvcat(str_doc_ini, tab_AO(21:ltab));



%  ligne 19 : Data manager
% ------------------------

	tab_AO   = fgetl(ftab_AO);
	ltab = length(tab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(19,:)
            h= warndlg(mess1,mess_arret);
            waitfor(h);
            break;
	end

	str_doc_ini= strvcat(str_doc_ini, tab_AO(21:ltab));



%  ligne 20 : Geographic Area
% ---------------------------

	tab_AO   = fgetl(ftab_AO);
	ltab = length(tab_AO);
	label_AO = tab_AO(1:20);
	if   label_AO ~= tab_label(20,:)
	    h= warndlg(mess1,mess_arret);
            waitfor(h);
            break;
	end

	str_doc_ini= strvcat(str_doc_ini, tab_AO(21:ltab));




     	 clear  tab_AO;

    	 fclose(ftab_AO);


     end     % fin du test sur le type de traitement (itr = 2, AO)



   else

% si problème de lecture de la table paramètres

     	    imes = 1;




end   % fin du test sur imess_choixtab


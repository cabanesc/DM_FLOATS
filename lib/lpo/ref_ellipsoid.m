function [a, b, N] = ref_ellipsoid(Name)
%
% 
% Return parameters of known ellipsoids
%
% WGS 84
% a = 6 378 137 m
% b= 6 356 752.3142 m
%
% IAG GRS 80
% L'ellipsoide du systeme International GRS 80 a ete adopté pour le système
% géodésique WGS 84 donc memes paramètres 
%
% Europe 50 Ellipsoide Hayford 1924 dit "International"
% a = 6 378 388
% b = 6 356 911.9
%
% Autre ellipsoide d'usage frequent Hayford 1909  : memes parametres que
% l'international
%
% Nouvelle Triangulation de la France  Ellipsoide Clarke 1880
% a = 6 378 249.1
% b = 6 356 514.9
%
% Dans les secteurs qui vous intéressent, l'ellipsoide Airy 1830, utilisé par
% la Grande Bretagne, est peut être en usage :
% a = 6 377 563.4 
% b = 6 356 256.9
%
%
% Pour memoire, f = 1 - b/a ou b = a(1-f)
%
% Source :
%           Eric MOUSSAT
%
% Banque de Bathymetrie & Geophysique Marine
% IFREMER/TMSI/IDM/SISMER BP70 29280 PLOUZANE FRANCE
% Bur.# (33) 02 98 22 42 05 Sec.# --- 45 41 Fax# --- 46 44  
% E-Mail : Eric.Moussat@ifremer.fr ou sismer@ifremer.fr                    
% 
% -------------------------------------
%

Ellipsoid_choice = ['WGS 84      '
			 		'IAG GRS 80  '
                    'Europa 50   '
                    'Hayford 1909'
                    'Clarke 1880 '
                    'Airy 1830   '];

Ellipsoid_choice_Mj_axis = [6378137.000
			    			6378137.000
			    			6378388.000
			    			6378388.000
			    			6378249.100
			    			6377563.400];

Ellipsoid_choice_mn_axis = [6356752.3142
			    			6356752.3142
                            6356911.9000
                            6356911.9000
                            6356514.9000
                            6356256.9000];                    

ix = strmatch(Name, deblank(Ellipsoid_choice));
if ~isempty(ix)
	a = Ellipsoid_choice_Mj_axis(ix);
	b = Ellipsoid_choice_mn_axis(ix);
	N = Ellipsoid_choice(ix,:);
else
	a = [];
	b = [];
	N = Ellipsoid_choice;
end	
return;

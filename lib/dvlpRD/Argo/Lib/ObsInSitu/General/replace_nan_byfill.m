function [Co]=replace_nan_by_fill(Co,fillvalName)
% -========================================================
%   USAGE : [Co]=replace_nan_by_fill(Co)
%           [Co]=replace_nan_by_fill(Co,fillvalName)
%
%   PURPOSE : short description of the function here
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
if isempty(strfind(Co.obj,'ObsInSitu'))
   error('replace_nan_byfill not define for this type of structure')
else 
    if nargin==1
       fillvalName='FillValue_';
    end
    champs = fieldnames(Co);    %champs={'psal','psalqc','psalad',....}
    Nbfields = length(champs);
    gotoend=0;
    if isfield(Co, 'fillisnan')
	if Co.fillisnan==0
	gotoend=1;
	end
    end
    
    if gotoend==0
	    for k=1:Nbfields            % boucle sur toutes les variables
	    
		oneChamp=champs{k};
		if isfield(Co.(oneChamp),'data')
		if isempty(Co.(oneChamp).data)==0    % test si la variable n'est pas vide
		
		if isnumeric(Co.(oneChamp).data)==1   % test si la variable est un tableau numerique
		
		    if isempty(Co.(oneChamp).(fillvalName))==0   % test si il y a une fillvalue pour cette variable
		    
			selec_fill = isnan( Co.(oneChamp).data);
			Co.(oneChamp).data(selec_fill) = Co.(oneChamp).(fillvalName);
			Co.fillisnan=0;
		    else
		    warning('Did not find a fillvalue attribute')
		    end
		end
		end
		end
	    end
    end
end
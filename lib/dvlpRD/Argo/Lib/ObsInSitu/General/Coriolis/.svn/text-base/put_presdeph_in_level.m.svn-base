function Pr = put_presdeph_in_level(pres, deph , Co, level)
% -========================================================
%   USAGE : Pr = put_presdeph_in_level_2(pres, deph , Co, level)
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
%   CALLED SUBROUTINES: find_level (interne)
% ========================================================
if isempty(strfind(Co.obj,'ObsInSitu/Coriolis'))
   error('put_presdeph_in_level not define for this type of structure')
else
    
    %%%% Tests de base
    Co = check_FirstDimArray_is(Co,'N_PROF');
    if isfield(Co,[pres])==0        % teste si le parametre pres existe
%	disp(['function: put_presdeph_in_level => does not find ', pres , ' field. Use ' deph ' only'])
    end
    if isfield(Co,[deph])==0
%	    disp(['function: put_presdeph_in_level => does not find ', deph , ' field. Use ' pres ' only'])
    end
    if isfield(Co,[pres])==0  &  isfield(Co,[deph])==0 
	error('did not find any pres or deph')
    end
    pres_in = pres;
    deph_in = deph;
    level_in = level;
    
    % Remplace par les profils "deph" tous les profils "pres" qui sont fillval sur tous les niveaux
    % => les profils donnant les niveaux sont dans Pr.level
    
    [Pr,isfillvalpres] = find_level ( pres_in, deph_in, Co, level_in);
    
    % isfillvalpres  est une variable logique qui repere les profils de pression qui sont fillval sur tous les niveaux
    %     =1 si tout le profil pres est == fillvalue 
    %     =0 si au moins une valeur differente de fillvalue pour le profil pres
    
    % Fait la meme chose pour toutes les variables associees (les flags qc et les variables ajustees)
    % si elles existent
    
    
    isfillvalpres_in=isfillvalpres; 
    
    %  pres_in=[pres 'qc'];
    %  deph_in=[deph 'qc'];
    %  level_in=[level 'qc'];
    %  
    %  if isfield(Pr, pres_in)
    %  [Pr] = find_level (pres_in, deph_in, Pr, level_in, isfillvalpres_in);
    %  end
    
    pres_in=[pres '_qc'];
    deph_in=[deph '_qc'];
    level_in=[level '_qc'];
    
    
    %keyboard
    %disp('toto')
    [Pr] = find_level (pres_in, deph_in, Pr, level_in, isfillvalpres_in);
   
    
    
    %  pres_in=[pres 'ad'];
    %  deph_in=[deph 'ad'];
    %  level_in=[level 'ad'];
    %  
    %  if  isfield(Pr, pres_in)
    %  [Pr] = find_level (pres_in, deph_in, Pr, level_in, isfillvalpres_in);
    %  end
    
    pres_in=[pres '_adjusted'];
    deph_in=[deph '_adjusted'];
    level_in=[level '_adjusted'];
    
    
    [Pr] = find_level (pres_in, deph_in, Pr, level_in, isfillvalpres_in);
  
    
    
    
    %  pres_in=[pres 'adqc'];
    %  deph_in=[deph 'adqc'];
    %  level_in=[level 'adqc'];
    %  
    %  if  isfield(Pr, pres_in)
    %  [Pr] = find_level (pres_in, deph_in, Pr, level_in, isfillvalpres_in);
    %  end
    
    pres_in=[pres '_adjusted_qc'];
    deph_in=[deph '_adjusted_qc'];
    level_in=[level '_adjusted_qc'];
    
    
    [Pr] = find_level (pres_in, deph_in, Pr, level_in, isfillvalpres_in);
   
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        SUBROUTINE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [Pr,varargout] = find_level ( pres, deph, Pri, level, isfillvalpres_in)
%keyboard
if nargin==4
   isfillvalpres_in=[];
end

Pr=Pri;
Pr.(level).data=[];

fillval='FillValue_';

%%%% On essaie de remplir level avec pres
%----------------------------------------
if isfield(Pri,[pres])        % teste si le parametre pres existe
    if isempty(Pri.(pres).data)==0; % si en plus il est plein, on remplie level 
	Pr.(level) = Pri.(pres);
	Pr.(level).isfrompres = logical(ones(size(Pri.(pres).data,1),1));
	isfillvalpres=logical(zeros(size(Pri.(pres).data,1),1));
	%disp([pres,' is filled']);
	
	if isempty(isfillvalpres_in)==0
		if isfield(Pr.(pres),fillval)==0; Pr.(pres).(fillval)=999.;end;
		Pr.(level).data(isfillvalpres_in==1,:) = repmat(Pr.(pres).(fillval),[sum(isfillvalpres_in==1),size(Pr.(pres).data,2)]);
	        Pr.(level).isfrompres=isfillvalpres_in==0;
	end
    end
end

if isfield(Pri,[deph])
    if isempty(Pri.(deph).data)==0; 
    
	if isempty(Pr.(level).data)==1
	   %disp([deph,' is filled']);
	%%%% Si level est toujours vide, on essaie de le remplir avec deph
	%-----------------------------------------------------------------
	    Pr.(level) = Pri.(deph);
	    Pr.(level).isfrompres = logical(zeros(size(Pri.(deph).data,1),1));
	    isfillvalpres=logical(ones(size(Pri.(deph).data,1),1));
	    if isempty(isfillvalpres_in)==0
		if isfield(Pr.(deph),fillval)==0; Pr.(deph).(fillval)=999.;end;
		Pr.(level).data(isfillvalpres_in==0,:) = repmat(Pr.(deph).(fillval),[sum(isfillvalpres_in==0),size(Pr.(deph).data,2)]);
	        Pr.(level).isfrompres=isfillvalpres_in==0;
	    end
	       
	else
	    %disp(['both ' deph ' and ' pres]);
	    %%%% Sinon  on regarde si les variables de niveau sont a la fois dans pres et deph (par exemple le cas pour fichiers _TE.nc)
	    % --------------------------------------------------------------------------------
	    
	    %%%% on essaie de remplacer par deph tous les profils qui ont pres = fillval sur tous les niveaux
	    % -----------------------------------------------------------------------------------------------
	    if isempty(isfillvalpres_in)==1
		isfillvalpres = logical(zeros(Pri.n_prof,1));
		if isempty (Pri.(pres).(fillval))==0 
		    isfillvalpres = sum((Pri.(pres).data == Pri.(pres).(fillval)),2)  ==  size(Pri.(pres).data,2);
		    % isfillvalpres =1 si tout le profil pres est == fillvalue 
		    % isfillvalpres =0 si au moins une valeur differente de fillvalue pour le profil pres
		end
	    else
		isfillvalpres=isfillvalpres_in;
	    end
	    
	    select_profils = (isfillvalpres==1);
	    select_profils;
	    Pr.(level).data(select_profils,:) = Pri.(deph).data(select_profils,:);
	    Pr.(level).isfrompres = (select_profils==0);
	    
	end
	
	if isempty(Pr.(level).data)==0 
	    Pr.(level).long_name=[Pri.(deph).long_name];  % on reste en unite de profondeur
	    Pr.(level).name=upper(level);
	    Pr.(level).units=[Pri.(deph).units];
        end
	
    end
end 


if isempty(Pr.(level).data)==1 
   Pr.(level).long_name=[];
   Pr.(level).name=[];
   Pr.(level).dim=[];
   Pr.(level).typ=[];
   Pr.(level).units=[];
   Pr.(level).(fillval)=[];

end
   

if nargout==2
    varargout{1}=isfillvalpres;
end
return
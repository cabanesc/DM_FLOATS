% -========================================================
%   USAGE : [OUT,..]=template_sub(IN,..)
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

function [Pr varargout] = adjust_param_and_qc (param,paramad,Co,newparam,MODE, select_profiles_in)

% param =le champ du parametre TR ex 'temp'
% paramad = le champ du parametre ajuste ex 'tempad'
% newparam = le champ qui contient le parametre TR remplace par DM quand c'est possible 

fillval='FillValue_';

%%% TEST DE BASE
Co = check_FirstDimArray_is(Co,'N_PROF');

rec_fillisnan=0;
if isfield(Co, 'fillisnan')
    rec_fillisnan=Co.fillisnan;
    if Co.fillisnan==1
    Co=replace_nan_byfill(Co,fillval);
    end
end

switch MODE
    case{'TR'}
    iadj=0;
    Pr.isadj='TR';
    case{'DM'}
    iadj=1;
    Pr.isadj='DM';
    MODE='D';
    case{'TRA'}
    iadj=1;
    Pr.isadj='TRA';
    MODE='A';
end

if isfield(Co,[param])==0
     disp(['function: adlust_param => does not find ', param , ' field'])
end
if isfield(Co,[paramad])==0
     disp(['function: adlust_param => does not find ', paramad , ' field'])
end



%%% INIT

select_profils=logical(zeros(size(Co.(param).data,1),1));

Pr=Co;



if isfield(Co,[param])        % teste si le parametre existe

    Pr.(newparam) = Co.(param);  % init par le champ TR
    Pr.(newparam).isadj = []; 

    if isempty(Co.(param).data)==0  % teste si le parametre n'est pas vide

	Pr.(newparam).isadj = logical(zeros(size(Co.(param).data,1),1));
	Pr.isadj='TR';
	
	if isfield(Co,[paramad])              % teste si le parametre ajuste existe
	
	    if iadj==1
	    
		if isempty(Co.(paramad).data)==0   % teste si le parametre ajuste n'est pas vide
		    % teste si le parametre ajuste a une fillvalue definie
		    
		    if isempty (Co.(paramad).(fillval))==0 
		    
			isfillvalprofile = sum((Co.(paramad).data == Co.(paramad).(fillval)),2)  ==  size(Co.(paramad).data,2);
			% isfillvalprofile =1 si tout le profil ajuste est == fillvalue 
			% isfillvalprofile =0 si au moins une valeur differente de fillvalue pour le profil ajuste 
		    end
		    
		    if nargin==5
		       select_profiles_in=logical(ones(size(Co.(param).data,1),1));
		    end
		    
		    select_profils = (isfillvalprofile==0)& select_profiles_in;
		    
		    if isempty(Co.dmode.data)==0
		    
			select_profils = select_profils &(Co.dmode.data==MODE); 
			
		        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if isempty(select_profils)==0
			
			    Pr.(newparam).data(select_profils,:) = Co.(paramad).data(select_profils,:);
			    
			    % les qc flags subissent le meme traitement
			    
			    if isfield(Co,[param,'_qc'])
			        Pr.([newparam '_qc']) = Co.([param 'qc']); % init par le champ TR
			        Pr.([newparam '_qc']).isadj = []; 
				if isempty(Co.([param,'_qc']).data)==0 
				   Pr.([newparam,'_qc']).isadj = select_profils;
				   
				    if  isfield(Co,[paramad,'_qc'])
					if isempty(Co.([paramad,'_qc']).data)==0 
						Pr.([newparam,'_qc']).data(select_profils,:) = Co.([paramad,'_qc']).data(select_profils,:);
						
					else
						Pr.([newparam,'_qc']).data(select_profils,:) = repmat(Co.([param,'_qc']).(fillval),[sum(select_profils),size(Co.([param,'_qc']).data,2)]);
					end
				   else
						Pr.([newparam,'_qc']).data(select_profils,:) = repmat(Co.([param,'_qc']).(fillval),[sum(select_profils),size(Co.([param,'_qc']).data,2)]);
				    end
				end
			    end
			    
			    %
			    Pr.([newparam]).isadj = select_profils;
			end
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		    end
		
		end
	    end
	    
        end
	
    end
end
if nargout==2
    varargout{1}=select_profils;
end


if rec_fillisnan==1
    Co=replace_fill_bynan(Co);
end

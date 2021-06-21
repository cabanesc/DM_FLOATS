function fprintf_to_latex_table(ficname)

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


nout=nargout;  % nombre de colonnes à lire (peut etre inferieur au nombre total de colonne)

if nargin==1
   D= ' ';
   select_ligne='';
end
if nargin==2
   select_ligne='';
end

except_ligne='#'; % ne lit pas les lignes debutant par '#'

icompt = 1; % determine la fin du fichier quand ==0
ilin=0;
nblin=0;
if exist(ficname,'file')==2
    fid=fopen(ficname);
    while icompt
	tline=fgetl(fid); % tline est de type char
	
	if ~ischar(tline) 
	    icompt=0;     % arrive a la fin du fichier
	else
	    if isspace(D)==0  % si jamais separateur de champ, remplace les separateurs accolés par deux separateur espacés
		
		tline= strrep(tline,[D D],[D ' ' D]);
		tline= strrep(tline,[D D],[D ' ' D]);
	    end
	    nblin=nblin+1;
	    
	    [part1,r]=strtok(tline,D); % coupe tline au niveau du premier espace (1ere col)  ' '=> part1 avant, r apres
	    if isempty(select_ligne)==1 % on lit toutes les lignes sans distinction 
		if strcmp(tline(1),except_ligne)==0 
		    ilin=ilin+1;            % compteur de ligne
		    num_ligne(ilin)=nblin;
		    col.l1{ilin}=part1; 
		    for nbcol=2:nout
			dynf=['l' num2str(nbcol)];
			[part1,r]=strtok(r,D);
			col.(dynf){ilin}=part1;
		    end
		end
	    else
		if findstr(part1,select_ligne)  % verifie que l'on est bien sur une ligne de modif de flag (qui commence par le nom du fichier netcdf)
		    if strcmp(tline(1),except_ligne)==0 
			ilin=ilin+1;         % compteur de ligne
			col.l1{ilin}=part1; 
			num_ligne(ilin)=nblin;
			for nbcol=2:nout
			    dynf=['l' num2str(nbcol)];
			    [part1,r]=strtok(r,D);
			    col.(dynf){ilin}=part1;
			end
		    end
		end
	    end
	end
    end
    fclose(fid);
else
    error(['le fichier ' ficname ' n''existe pas'])
end
for nbcol=1:nout
    dynf=['l' num2str(nbcol)];
    varargout{nbcol,:} = {col.(dynf){:}};
end
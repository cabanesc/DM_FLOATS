function [date_juld]=extract_date_file_name(filename)
% -========================================================
%   USAGE : [date_juld]=extract_date_file_name(filename)
%   PURPOSE : extrait la date du nom d'un fichier Coriolis
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

% specifique aux fichiers CORIOLIS avec date de la forme 'yyyymmdd': A AMELIORER

% extrait la date du nom du fichier
islash = findstr(filename,'/');

fic = filename(islash(end)+1:end);

% trouve l'annee et donc le repertoire

% trouve 8 chiffres consecutifs
A=  isstrprop(fic, 'digit');
B=(A(1:end-7)==1&A(2:end-6)==1&A(3:end-5)==1&A(4:end-4)==1&A(5:end-3)==1&A(6:end-2)==1&A(7:end-1)==1&A(8:end)==1);
kb=find(B==1);
if length(kb)~=1
    error('Verifier le nom des fichiers : ils ne doivent contenir qu une (seule) date yyyymmdd')
else
    annee=fic(kb:kb+3) ;
    mois=fic(kb+4:kb+5);
    jour=fic(kb+6:kb+7) ;
end  
date_juld = datenum([annee '-' mois '-' jour], 'yyyy-mm-dd')-datenum('1950-01-01','yyyy-mm-dd');
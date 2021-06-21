function fprintf_to_latex_table_3(fw,ind,tabular,lignetowrite)

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

if nargin==3
   lignetowrite='';
elseif nargin<3
   error(' fprintf_to_latex_table : Not enough imput argument')
end

if ind==0 % entete du fichier latex
    fprintf(fw,'%s\n', [' \documentclass[12pt]{article}']);
    fprintf(fw,'%s\n', ['\usepackage{longtable}']);
    fprintf(fw,'%s\n', ['\usepackage{array}']);
    fprintf(fw,'%s\n', ['\begin{document}']);
    fprintf(fw,'%s\n', ['\newcolumntype{x}[1]{%']);
    fprintf(fw,'%s\n', ['>{\raggedright}p{#1}}%']);
    fprintf(fw,'%s\n', ['\begin{center}']);
    fprintf(fw,'%s\n', ['\begin{longtable}{' tabular '}']);
    fprintf(fw,'%s\n', ['\hline']);
    
    
    if isempty(lignetowrite)==0
	lignetowrite=strrep(lignetowrite,',',' & ');
	lignetowrite=strrep(lignetowrite,'\_','_');
	lignetowrite=strrep(lignetowrite,'_','\_');
	lignetowrite=strtrim([lignetowrite '\tabularnewline']);
	fprintf(fw,'%s\n',lignetowrite);
	fprintf(fw,'%s\n','\hline');
    end
    fprintf(fw,'%s\n', ['\hline']);
    fprintf(fw,'%s\n', ['\endfirsthead']);
    fprintf(fw,'%s\n', ['\multicolumn{4}{c}%']);
    fprintf(fw,'%s\n', ['{\tablename\ \thetable\ -- \textit{Continued from previous page}} \\']);
    fprintf(fw,'%s\n', ['\hline']);
    if isempty(lignetowrite)==0
	fprintf(fw,'%s\n',lignetowrite);
	fprintf(fw,'%s\n', ['\hline']);
    end
    fprintf(fw,'%s\n', ['\endhead']);
    fprintf(fw,'%s\n', ['\hline \multicolumn{4}{r}{\textit{Continued on next page}} \\']);
    fprintf(fw,'%s\n', ['\endfoot']);
    fprintf(fw,'%s\n', ['\hline']);
    fprintf(fw,'%s\n', ['\endlastfoot']);
end

if ind==1   % corp du tableau
   if isempty(lignetowrite)==0
       lignetowrite=strrep(lignetowrite,',',' & ');
	lignetowrite=strrep(lignetowrite,'\_','_');
	lignetowrite=strrep(lignetowrite,'_','\_');
	lignetowrite=[lignetowrite '\tabularnewline'];
	fprintf(fw,'%s\n',lignetowrite);
	%fprintf(fw,'%s\n','\hline');
    
   end
end

if ind==2   % corp du tableau
    if isempty(lignetowrite)==0
	lignetowrite=strrep(lignetowrite,',',' & ');
	lignetowrite=strrep(lignetowrite,'\_','_');
        lignetowrite=strrep(lignetowrite,'_','\_'); 
	lignetowrite=[lignetowrite '\tabularnewline'];
	fprintf(fw,'%s\n',lignetowrite);
	%fprintf(fw,'%s\n','\hline');
    end
    fprintf(fw,'%s\n',['\end{longtable}']);
    fprintf(fw,'%s\n',['\end{center}']);
    fprintf(fw,'%s\n',['\end{document}']);
    
end


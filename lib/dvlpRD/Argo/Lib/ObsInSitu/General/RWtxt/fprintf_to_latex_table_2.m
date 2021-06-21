function fprintf_to_latex_table(fw,ind,tabular,lignetowrite)

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

    fprintf(fw,'%s\n',['\documentclass[a4paper,11pt]{article}']);
    
    
    fprintf(fw,'%s\n',['\usepackage{tabularx}']);
    fprintf(fw,'%s\n',['\usepackage{longtable}']);
    fprintf(fw,'%s\n',['\usepackage[dvips,final]{graphicx}']);
    fprintf(fw,'%s\n',['\usepackage{lscape}']);
    fprintf(fw,'%s\n',['\usepackage{color}']);
    fprintf(fw,'%s\n',['\usepackage{colortbl}']);
    fprintf(fw,'%s\n',['\usepackage{hyperref}']);
    fprintf(fw,'%s\n',['\hypersetup{urlcolor=blue,linkcolor=blue,citecolor=black,colorlinks=true}']);
    
    fprintf(fw,'%s\n',['\begin{document}']);
    
    fprintf(fw,'%s\n',['\begin{table}[h!]']);
    fprintf(fw,'%s\n',['{%']);
    fprintf(fw,'%s\n',['\renewcommand{\arraystretch}{1.5}']);
    fprintf(fw,'%s\n',['\newcommand{\mc}[3]{\multicolumn{#1}{#2}{#3}}']);
    fprintf(fw,'%s\n',['\newcolumntype{M}[1]{>{\raggedright}m{#1}}']);

    fprintf(fw,'%s\n',['\definecolor{tcA}{rgb}{0.917647,0.917647,0.917647}']);
    fprintf(fw,'%s\n',['\begin{center}']);
    fprintf(fw,'%s\n',['\small']); 
    fprintf(fw,'%s\n',['\begin{longtable}{' tabular '}']);
    fprintf(fw,'%s\n',['\rowcolor{tcA}']);
    fprintf(fw,'%s\n',['\hline']);
    if isempty(lignetowrite)==0
	lignetowrite=strrep(lignetowrite,',',' & ');
	lignetowrite=strrep(lignetowrite,'\_','_');
	lignetowrite=strrep(lignetowrite,'_','\_');
	lignetowrite=strtrim([lignetowrite '\\']);
	fprintf(fw,'%s\n',lignetowrite);
	fprintf(fw,'%s\n','\hline');
    end
end

if ind==1   % corp du tableau
   if isempty(lignetowrite)==0
       lignetowrite=strrep(lignetowrite,',',' & ');
	lignetowrite=strrep(lignetowrite,'\_','_');
	lignetowrite=strrep(lignetowrite,'_','\_');
	lignetowrite=[lignetowrite '\\'];
	fprintf(fw,'%s\n',lignetowrite);
	fprintf(fw,'%s\n','\hline');
    
   end
end

if ind==2   % corp du tableau
    if isempty(lignetowrite)==0
	lignetowrite=strrep(lignetowrite,',',' & ');
	lignetowrite=strrep(lignetowrite,'\_','_');
        lignetowrite=strrep(lignetowrite,'_','\_'); 
	lignetowrite=[lignetowrite '\\'];
	fprintf(fw,'%s\n',lignetowrite);
	fprintf(fw,'%s\n','\hline');
    end
    fprintf(fw,'%s\n',['\end{longtable}']);
    fprintf(fw,'%s\n',[' \normalsize']);
    fprintf(fw,'%s\n',['\label{tab4}']);
    fprintf(fw,'%s\n',['\end{center}']);
    fprintf(fw,'%s\n',[' }%']);
    fprintf(fw,'%s\n',['\end{table}']);
    fprintf(fw,'%s\n',['\end{document}']);
    
end


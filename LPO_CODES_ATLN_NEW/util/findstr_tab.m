function [isfound]=findstr_tab(tab,car)
% -========================================================
%   USAGE : [isfound]=findstr_tab(TAB,car)
%   PURPOSE : trouve une chaine de caractère dans un tableau de caractere
%   ou cell array
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

if iscell(tab)==0 
    % transforme en cell
    tabcell=cellstr(tab);
else
    tabcell=tab;
end
% 2 on cherche la chaine de caratere
ischarcell = strfind(tabcell,car );

% on cherche ou la cell ischarcell n'est aps vide
isfound=~cellfun('isempty',ischarcell);

    

function [Coi]=extract_profile(Co,DIMNAME,iprofiles)
% -========================================================
%   USAGE : [Coi]=extract_profile(Co,DIMNAME,iprofiles)
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
    error('extract_profile not define for this type of structure')
else
    INITFIRSTDIM=[];
    if isfield(Co,'firstdimname')
        INITFIRSTDIM=Co.firstdimname;
    end
    Co = check_FirstDimArray_is(Co,DIMNAME);
    Coi=Co;

    champs = fieldnames(Co);    %champs={'psal','psalqc','psalad',....}
    Nbfields = length(champs);

    for k=1:Nbfields            % boucle sur toutes les variables
        oneChamp=champs{k};
        if isfield(Co.(oneChamp),'data')
            if isempty(Co.(oneChamp).data)==0
                isthedim=strcmp(Co.(oneChamp).dim,DIMNAME);
                if sum(isthedim)==1
                    nbdim = length(size(Coi.(oneChamp).data));
                    ap='';
                    if nbdim>1
                        ap=repmat(',:',[1,nbdim-1]);
                    end
                    expre=['Coi.(oneChamp).data=Co.(oneChamp).data(iprofiles' ap ');'];
                    eval(expre);
                    %Coi.(oneChamp).data=Co.(oneChamp).data(iprofiles,:);
                end
            end
        end
    end

    Coi = check_FirstDimArray_is(Coi,DIMNAME);
    if isempty(INITFIRSTDIM)==0
        Coi = check_FirstDimArray_is(Coi,INITFIRSTDIM);
    end
end

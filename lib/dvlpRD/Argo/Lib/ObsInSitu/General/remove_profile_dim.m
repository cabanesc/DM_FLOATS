function [Coi,Dimi]=remove_profile(Co,Dim,DIMNAME,iprofiles)
% -========================================================
%   USAGE : [Coi]=remove_profile(Co,DIMNAME,iprofiles)
%   EXAMPLE [Coi]=remove_profile(Co,'N_PROF',[4,8,12])
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
   error('remove_profile not define for this type of structure')
else
    if islogical(iprofiles)==1|sum((iprofiles==0|iprofiles==1))==length(iprofiles)
       iprofiles=find(iprofiles==1);
    end
    if isfield(Co,'firstdimname')
    INITFIRSTDIM=Co.firstdimname;
    else
    INITFIRSTDIM=[];
    end
    
    Co = check_FirstDimArray_is(Co,DIMNAME);
    Coi=Co;
    Dimi=Dim;

    champs = fieldnames(Co);    %champs={'psal','psalqc','psalad',....}
    Nbfields = length(champs);
    
    for k=1:Nbfields            % boucle sur toutes les variables
	oneChamp=champs{k};
	if isfield(Co.(oneChamp),'data')
	if isempty(Co.(oneChamp).data)==0
	    isthedim=strcmp(Co.(oneChamp).dim,DIMNAME);
	    if sum(isthedim)==1
		    allprofiles=[1:size(Coi.(oneChamp).data,1)];
		    keepprofiles= setdiff(allprofiles,iprofiles);
		    nbdim = length(size(Coi.(oneChamp).data));
		    ap='';
		    if nbdim>1
			ap=repmat(',:',[1,nbdim-1]);
		    end
		    expre=['Coi.(oneChamp).data=Co.(oneChamp).data(keepprofiles' ap ');'];
		    eval(expre);
	    end
	end
	end
    end
    dimnamei=fieldnames(Dimi);
   % keyboard
    for k=1:length(dimnamei)
           if strcmp(dimnamei{k},lower(DIMNAME))==1 % dimension selon laquelle les champs sont extraits
           Dimi.(dimnamei{k}).dimlength=Dim.(dimnamei{k}).dimlength-length(iprofiles);
           end
    end
    
    Coi = check_FirstDimArray_is(Coi,DIMNAME);
    if isempty(INITFIRSTDIM)==0
	Coi = check_FirstDimArray_is(Coi,INITFIRSTDIM);
    end
end
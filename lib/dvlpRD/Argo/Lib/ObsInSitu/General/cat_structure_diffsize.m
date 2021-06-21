function [Coi]=cat_structure(Co1,Co2,DIMNAME)
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
% TODO: A MODIFIER pour tenir compte des cas ou les structures n'ont pas la meme taille sur une dimension(voir replace profile)
if isfield(Co1,'firstdimname')
    INITFIRSTDIM1=Co1.firstdimname;
else
    INITFIRSTDIM1=[];
end
if isfield(Co2,'firstdimname')
    INITFIRSTDIM2=Co2.firstdimname;
else
    INITFIRSTDIM2=[];
end
    
if isempty(Co1)==1
   Coi=Co2;
else
    Co1 = check_FirstDimArray_is(Co1,DIMNAME);
    Co2 = check_FirstDimArray_is(Co2,DIMNAME);
    
    ndim1 = Co1.(lower(DIMNAME));
    ndim2 = Co2.(lower(DIMNAME));
    
    Coi=Co1;
    
    champs1 = fieldnames(Co1);    %champs={'psal','psalqc','psalad',....}
    champs2 = fieldnames(Co2);
    
    
    Nbfields1 = length(champs1);
    Nbfields2 = length(champs2);
    
    for k=1:Nbfields1            % boucle sur toutes les variables
	oneChamp=champs1{k};
	if isfield(Co2,oneChamp)
	if isfield(Co1.(oneChamp),'data')&isfield(Co2.(oneChamp),'data')
	if isempty(Co1.(oneChamp).data)==0&isempty(Co2.(oneChamp).data)==0
	
	    isthedim1=strcmp(Co1.(oneChamp).dim,DIMNAME);
	    isthedim2=strcmp(Co2.(oneChamp).dim,DIMNAME);
	    if sum(isthedim1)==1&sum(isthedim2)==1
	        otherdim1=(isthedim1==0);
		otherdim2=(isthedim2==0);
		fullsize1 = size(Co1.(oneChamp).data);
		size1 = fullsize1(otherdim1);
		fullsize2 = size(Co2.(oneChamp).data);
		size2 = fullsize2(otherdim2);
		maxsize=max(size1,size2);
		nbdim1=length(size1);
		nbdim2=length(size2);
		if isequal(size1,maxsize)==0
		    ap='';
		    for m=1:nbdim1
			ap=[ap ',1:size1(' num2str(m) ')'];
		    end
		    themax=[fullsize1(1),maxsize];
		    if isfield(Co1.(oneChamp),'FillValue_')
		    Co.(oneChamp).data=repmat(Co1.(oneChamp).FillValue_,themax);
		    expre=['Co.(oneChamp).data(:' ap ')=Co1.(oneChamp).data;'];
		    eval(expre)
		    else
		    oneChamp
		    error('Did not find FillValue for this field')
		    end
		    Co1.(oneChamp).data=Co.(oneChamp).data;
		    clear Co
		end
		if isequal(size2,maxsize)==0
		    ap='';
		    for m=1:nbdim2
			ap=[ap, ',1:size2(' num2str(m) ')'];
		    end
		    themax=[fullsize2(1),maxsize];
		    tempo = Co2.(oneChamp).data;
		    if isfield(Co2.(oneChamp),'FillValue_')
		    Co2.(oneChamp).data = repmat(Co2.(oneChamp).FillValue_,themax);
		    expre=['Co2.(oneChamp).data(:' ap ')=tempo;'];
		    eval(expre)
		    clear tempo
		    else
		    oneChamp
		    error('Did not find FillValue for this field')
		    end
		end
		Coi.(oneChamp)=Co1.(oneChamp);
		% allprofiles=[1:size(Coi.(oneChamp).data,1)];
		% keepprofiles= setdiff(allprofiles,iprofiles);
		% verifie que tous les champs "attributs" sont les memes
		attributs1 = rmfield(Co1.(oneChamp),'data');
		attributs2 = rmfield(Co2.(oneChamp),'data');
		if isequal(attributs1,attributs2)==0
		    error('The 2 structures have not the same attributes');
		else
		    nbdim = length(size(Coi.(oneChamp).data));
		    ap='';
		    if nbdim>1
		    ap=repmat(',:',[1,nbdim-1]);
		    end
		    expre=['Coi.(oneChamp).data(ndim1+1:ndim1+ndim2' ap ') = Co2.(oneChamp).data(:' ap ');'];
		    eval(expre);
		end
	    end
	end
	end
	end
    end

end

Coi = check_FirstDimArray_is(Coi,DIMNAME);
if isempty(INITFIRSTDIM1)==0
    Coi = check_FirstDimArray_is(Coi,INITFIRSTDIM1);
end

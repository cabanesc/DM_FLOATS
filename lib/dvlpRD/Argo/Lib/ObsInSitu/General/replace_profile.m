function [Co]=replace_profile(Coi,Cor,DIMNAME,iprofiles,rprofiles)
% -========================================================
%   USAGE : [Co]=replace_profile(Coi,Cor,DIMNAME,iprofiles,rprofiles)
%   EXAMPLE : [Co]=replace_profile(Coi,Cor,'N_PROF',[3,4,5],[6,7,8])
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
if isempty(strfind(Coi.obj,'ObsInSitu'))
   error('replace_profile not define for this type of structure')
else

    if length(iprofiles)~=length(rprofiles)
    error('iprofiles and rprofiles should be vector of same length')
    end

    INITFIRSTDIM=[];
    if isfield(Coi,'firstdimname')
	INITFIRSTDIM=Coi.firstdimname;
    end
    
    Cor = check_FirstDimArray_is(Cor,DIMNAME);
    Coi = check_FirstDimArray_is(Coi,DIMNAME);
    Co=Coi;
    
    champs = fieldnames(Coi);    %champs={'psal','psalqc','psalad',....}
    Nbfields = length(champs);
    
    for k=1:Nbfields            % boucle sur toutes les variables
	oneChamp=champs{k};
	if isfield(Cor,oneChamp)
	if isfield(Coi.(oneChamp),'data') & isfield(Cor.(oneChamp),'data')
	if isempty(Coi.(oneChamp).data)==0 & isempty(Cor.(oneChamp).data)==0
	if sum(strcmp(Coi.(oneChamp).dim ,Cor.(oneChamp).dim))==length(Coi.(oneChamp).dim)  % les champs ont les memes dim names
	    isthedim_i=strcmp(Coi.(oneChamp).dim,DIMNAME);
	    isthedim_r=strcmp(Cor.(oneChamp).dim,DIMNAME);
	    
	    if sum(isthedim_i)==1&sum(isthedim_r)==1
		otherdim_i=(isthedim_i==0);
		otherdim_r=(isthedim_r==0);
		fullsizei = size(Coi.(oneChamp).data);
		sizei = fullsizei(otherdim_i);
		fullsizer = size(Cor.(oneChamp).data);
		sizer = fullsizer(otherdim_r);
		maxsize=max(sizei,sizer);
		nbdim_i=length(sizei);
		nbdim_r=length(sizer);
		
		
		if isequal(sizei,maxsize)==0
		    ap='';
		    for m=1:nbdim_i
			ap=[ap ',1:sizei(' num2str(m) ')'];
		    end
		    themax=[fullsizei(1),maxsize];
		    if isfield(Coi.(oneChamp),'FillValue_')
		    Co.(oneChamp).data=repmat(Coi.(oneChamp).FillValue_,themax);
		    expre=['Co.(oneChamp).data(:' ap ')=Coi.(oneChamp).data;'];
		    eval(expre)
		    else
		    oneChamp
		    error('Did not find FillValue for this field')
		    end
		end
		
		if isequal(sizer,maxsize)==0
		    ap='';
		    for m=1:nbdim_r
			ap=[ap, ',1:sizer(' num2str(m) ')'];
		    end
		    themax=[fullsizer(1),maxsize];
		    tempo = Cor.(oneChamp).data;
		    if isfield(Cor.(oneChamp),'FillValue_')
		    Cor.(oneChamp).data = repmat(Cor.(oneChamp).FillValue_,themax);
		    expre=['Cor.(oneChamp).data(:' ap ')=tempo;'];
		    eval(expre)
		    clear tempo
		    else
		    oneChamp
		    error('Did not find FillValue for this field')
		    end
		end
		
		nbdim = length(size(Co.(oneChamp).data));
		ap='';
		if nbdim>1
		    ap=repmat(',:',[1,nbdim-1]);
		end
		expre=['Co.(oneChamp).data(iprofiles' ap ') = Cor.(oneChamp).data(rprofiles' ap ');'];
		eval(expre)
		
	    end
	end
	end
	end
	end
    end
    Co = check_FirstDimArray_is(Co,DIMNAME);
    if isempty(INITFIRSTDIM)==0
	Co = check_FirstDimArray_is(Co,INITFIRSTDIM);
    end
end
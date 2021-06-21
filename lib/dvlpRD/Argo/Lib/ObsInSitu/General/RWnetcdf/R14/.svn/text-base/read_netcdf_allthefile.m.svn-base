function [S,Dim,Globatt]=read_netcdf_allthefile(Ficname,VarIN,verbose)
% -========================================================
%   PURPOSE : read variables and attributes from netcdf files and put into a structure
%   USAGE : [S,Dim,Globatt]=read_netcdf_allthefile(Ficname)  % read all the variables
%           [S,Dim,Globatt]=read_netcdf_allthefile(Ficname,VarIN) % read only VarIN variables
%           [S,Dim,Globatt]=read_netcdf_allthefile(Ficname,1)  % print to debug
%
% -----------------------------------
%   INPUT :
%     Ficname   (string)  name of the file              ex: 'filetoread.nc'
%
%   OPTIONNAL INPUT :
%    VarIN  (structure)  name of variables to read      ex: To read only the variables PRES and TEMP in the file                          %                                                        set:                       
%                                                        VarIN.pres.name = 'PRES';
%                                                        VarIN.temp.name ='TEMP';
%                                                        rem: the fields for the structure VarIN are lower(varname)
%                        default: read all variables
%
%    verbose (scalar)    0 (default),  1 print to screen
% -----------------------------------
%   OUTPUT :
%    S (structure)
%        The structure S contains all info found in the netcdf file for the variable (including attributes), plus:
%        S.(lower(varname)).name   (string)   variable name       S.temperature.name = 'TEMPERATURE'                               
%        S.(lower(varname)).dim    (cell array of string)         S.temperature.dim  = {'N_PROF','N_LEVEL'}         
%        S.(lower(varname)).data        (array)                   S.temperature.data =  n_prof x n_level value of temperature
%
%
%     Dim   (structure)  
%        The structure Dim contains dimension info found in the netcdf file
%         Dim.n_prof.name = 'N_PROF'   (name of the dimension)
%         Dim.n_prof.dimlength  = 88   (length of the dimension)
% 
%     Globatt  (structure)  
%        The structure Globatt contains all the global attributes found in the netcdf file
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%   CALLED SUBROUTINES: none
% ========================================================

if nargin==1
    verbose=0;
    VarIN=[];
end

if nargin==2
    if isstruct(VarIN)==1
        verbose=0;
    else
        verbose=VarIN;
        VarIN=[];
    end
end

% Possible names for the fillvalue attribute
poss_fillval_name={'FillValue','FillValue_','_FillValue','_fillvalue','fillvalue','missing_value'};

% 1. Open the netcdf file
nc=netcdf(Ficname,'nowrite');

% test sur l'ouverture du fichier
if(isempty(nc)==1)
    error('===== > get_netcdf_files: Check the netcdf file name (it does not exist)')
end


thevars=var(nc);
Nbfields=size(thevars,2);

thedims = dim(nc);
Nbdims = size(thedims,2);

theglobatt = att(nc);
Nbglob = size(theglobatt,2);

% Donne la dimension 'unlimitted' ou record dimension
theRecdim=recdim(nc);


% recupere toutes les dimensions dans la structure Dim
Dim=[];
for kd=1:Nbdims
    namedim=name(thedims{kd});
    onedim=lower(namedim);
    Dim.(onedim).name=namedim;
    Dim.(onedim).dimlength=sum(size(thedims{kd}))-1;
end

% recupere toutes les global attributes dans la structure Globatt
Globatt=[];
for kd=1:Nbglob
    nameatt=name(theglobatt{kd});
    Globatt.(nameatt).att=theglobatt{kd}(:);
    Globatt.(nameatt).name=nameatt;
end

% Determine quelle sont les variables à lire
kread=[1:Nbfields];
for k=kread                  % trouve tous les noms des variables
    allvar{k}=deblank(name(thevars{k}));
end
inb=0;
if isempty(VarIN)==0
    tempo=fieldnames(VarIN);
    for nb=1:length(tempo)
        if isfield(VarIN.(tempo{nb}),'name')
            inb=inb+1;
            varTOretrieve{inb}=VarIN.(tempo{nb}).name;
        end
    end
    [com,ia,ib]=intersect(allvar,varTOretrieve);
    kread=kread(ia);
end

ndim=1;
% 2. Fill the structure S for each variable
for k=kread                  % boucle sur chaque variable a trouver ds le fichier

    VarName=deblank(name(thevars{k}));
    oneChamp=lower(VarName);


    clear l
    l=nc{VarName};                % on recupere la variable et tous ces attributs
   

    S.(oneChamp).name=VarName;
    S.(oneChamp).dim=[];          % => initialisation :tableau vide
    S.(oneChamp).data=[];
   % S.(oneChamp).units=[];

    vecdim=[];

    % Look for  attributes:
    % recherche tous les attributs de la variable
    attributs=att(l);                       % (att)
    clear name_att
    for iat=1:length(attributs)
        name_att{iat}=name(attributs{iat}); % recupere leur nom % outils nctools (name)
        % initialisation de la structure S, avec tous les attributs (sauf le fillvalue qui est remplit plus tard)
        isthefil=0;
        for ipos=1:length(poss_fillval_name)
            if strcmp(name_att{iat},poss_fillval_name{ipos});
                isthefil=1;
            end
        end
        if isthefil==0
            expre=['S.(oneChamp).' name_att{iat} '=l.' name_att{iat} '(:);'];
            eval(expre);
        else
            S.(oneChamp).FillValue_=[];
            % la recuperation du fillvalue se fait apres
        end
    end

    % recupere les dimensions
    dimen=dim(l);                    
    S.(oneChamp).dim='';

    % met le nom des dimensions dans une cellule de charactere pour y acceder facilement
    charToEval='{';                      % characters string
    for idim=1:length(dimen)             % work on each dimension
        thedimname=name(dimen{idim});    % name of the dimension  
        alldim{ndim}=thedimname;ndim=ndim+1;
        charToEval=[charToEval,'''', thedimname,'''',','];
        vecdim(idim)=sum(size(dimen{idim}))-1; % vecteur contenant les dimensions de la variable (il vaut mieux la calculer comme ca que par length!)
    end                                        % rem: length([0x1])= 1 alors qu'ici on voudrait 0

    charToEval=[charToEval,'}'];
    charToEval=strrep(charToEval,',}','}');

    S.(oneChamp).dim=eval(charToEval);
    % recupere le type  ex 'float' 'string8' ...
    S.(oneChamp).type=datatype(l); var_type=datatype(l);  
    % recupere les donnees
    S.(oneChamp).data=l(:);

    % recupere le FIllValue
    % Note: matlab cannot put the "_FillValue" attribute as a field in a structure (because of the "_" )
    % however matlab managed to read this attribute and put it in the "FillValue_" field

    % Look for missing value attribute:
    % recherche tous les attributs de la variable
    attributs=att(l);                       % (att)
    clear name_att
    missval=[];
    if isempty(attributs)==0
        for iat=1:length(attributs)
            name_att{iat}=name(attributs{iat}); % recupere leur nom (name)
        end
        ipos=1;
        while isempty(missval) & ipos<length(poss_fillval_name)
            isthefillatt=strcmp( name_att,poss_fillval_name{ipos});
            index=find(isthefillatt==1);
            if isempty(index)==0
                missval=l.(name(attributs{index}))(:);
            end
            ipos=ipos+1;
        end
    end

    if isempty(missval)==0
        S.(oneChamp).FillValue_=missval;
    end

    % Reshape the S.(oneChamp).data to coincide with vecdim
    if sum(vecdim==0)==0
        if sum(vecdim==1)>0 % presence de dimension singleton
            vecdim=[vecdim,1]; % evite un bug si le champ est scalaire
            S.(oneChamp).data=reshape(S.(oneChamp).data,vecdim);
        end
    end

    if verbose==1
        disp(['............',oneChamp])
    end
    if isfield(S.(oneChamp),'units')==0
    S.(oneChamp).units=[];
    end
end

% Clean not used dimensions for variables that are not retrieved

namedim=fieldnames(Dim);
rmdim=setdiff(namedim,lower(alldim));
for jdim=1:length(rmdim)
    Dim=rmfield(Dim,rmdim{jdim});
end

if isempty(theRecdim)==0
    if isempty(rmdim)==0
        if strcmp(lower(name(theRecdim)),rmdim)==0
            S.recdim=name(theRecdim);
        end
    else
        S.recdim=name(theRecdim);
    end
end

S.obj='ObsInSitu';
if isfield(VarIN,'obj')
    S.obj=[S.obj '/' VarIN.obj];
end
close(nc)

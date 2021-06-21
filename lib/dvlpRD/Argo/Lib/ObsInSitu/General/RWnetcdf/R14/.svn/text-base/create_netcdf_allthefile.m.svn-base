%_____________________________________________________________________________________________________
%   USAGE create_netcdf_allthefile(DIMS,VARS,ficout,CONFIG,verbose)
%%
%   PURPOSE: Create a netcdf file containing variables in VARI
%
%-----------------------------
%  INPOUT :
%
%   DIMS (structure)  - 
%       with required fields for each dim : .name        (string)           Dim.n_prof.name = 'N_PROF'     name of the dimension
%                                           .dimlength   (scalar)           Dim.n_prof.dimlength = 88      length of the dimension
%  VARS (structure)
%        with required fields for each variable:
%                               .name        (string)                       S.temp.name = 'TEMPERATURE' name of the variable
%                               .dim         (cell array of string)         S.temp.dim  = {'N_PROF','N_LEVEL'} 
%                               .data        (array)                        S.temp.data =  n_prof x n_level value of temperature
%
%        ficout (string) -  output file name-
%
%        CONFIG (structure)   OPTIONNAL contains global attribute you want to put in the netcdf file
%        with required fields : .name        (string)                        CONFIG.version.name='version';
%                             : .att         valeur de l'attribut global     CONFIG.version.att='1.1';
%-----------------------------
%  OUTPUT :
%   none
%-----------------------------
%  HISTORY :created mai 3007 ccabanes
%           revision oct2008, fev2009 ccabanes
%  CALLED SUBROUTINES:
%   define_var (internal)
%   store_var (internal)
%  ________________________________________________________________________________________________
function create_netcdf_allthefile(DIMS,VARS,ficout,CONFIG,verbose)
if nargin<=3
    CONFIG=[];
    verbose=0;
end

if nargin==4
    if isstruct(CONFIG)==1
        verbose=0;
    else
        verbose=CONFIG;
        CONFIG=[];
    end
end
% trouve la dimension unlimitted
REC=[];
REC_length=[];
if isstruct(VARS)==1
    if isfield(VARS,'recdim')
        if isempty(VARS.recdim)==0
            REC=VARS.recdim;
            REC_length=DIMS.(lower(REC)).dimlength;
        end
    end
end

Dim=DIMS;
% transforme en cellules
if isstruct(VARS)==1
    VARS=struct2cell(VARS);
end

if isstruct(DIMS)==1
    DIMS=struct2cell(DIMS);
end

NBVAR1=length(DIMS);
NBVAR2=length(VARS);
NBVAR=NBVAR1+NBVAR2;

for k=1:NBVAR1
    VARI{k}=DIMS{k};
end
if verbose==1
    VARI
end
for k=1:NBVAR2
    ik=k+NBVAR1;
    VARI{ik}=VARS{k};
end

clear VARS DIMS

%% Open the file:
%if verbose==1
%disp('=====================================================')
disp(['writting ' num2str(NBVAR2), ' variable(s)  in the file....'])
disp(ficout)
%end
f = netcdf(ficout, 'clobber');% Create netcdf file, as much variables as you want but same grid info!

% _________________________________________________________________________________________________________________________
%% Global attributes:

% Check in CONFIG if the fields exist

if isempty(CONFIG)==0
    if isstruct(CONFIG)==1
        CONFIG=struct2cell(CONFIG);
    end
    NBATT=length(CONFIG);
    for k=1:NBATT
        expre=['f.' CONFIG{k}.name ' = ncchar( CONFIG{k}.att);'];
        if verbose==1
            disp(['Global attribute : ' expre])
        end
        eval(expre)
    end
end


memo=whos('VARI'); % look for the memory space needed to store VARI -> memo.bytes

%% Dimensions et variables:
for k=1:NBVAR
   
    VARI{k};
    % init
    ADD_OFFSET=[];
    SCALE_FACTOR=[];
    missval=[]; %
    clear expr
    expr{1}=[]; expr{2}=[];
    expr{3}=[];expr{4}=[];
    expr{5}=[]; expr{6}=[];
    expr{7}=[];expr{8}=[];



    if isfield(VARI{k},'name')==0 & k<=NBVAR1
        disp(VARI{k})
        error ('should give a name  in Var.name to all dimensions')
    end

    if verbose==1
        if isfield(VARI{k},'name')
            disp(' ')
            disp(['............',VARI{k}.name])
        end
    end

    % Collecte les infos necessaires pour la suite si c'est une variable
    % ------------------------------------------------------------------
    if isfield(VARI{k},'data')==1
        % si c'est une variable (contient un champ .data)

        if isfield(VARI{k},'dim')==0
             disp(VARI{k})
            error ('should give a dim  in Var.dim to all variables')
        end

        % the dimension of the variable
        extract_dim=VARI{k}.dim;

        % screen output
        if isfield(CONFIG,'missval')
            missval=CONFIG.missval;
        end
        if isfield(VARI{k},'FillValue_')
            missval=VARI{k}.FillValue_;
        end
        if isfield(VARI{k},'fillval')
            missval=VARI{k}.fillval;
        end

        if isfloat(VARI{k}.data)|isinteger(VARI{k}.data) % single or double

            % check for infinite value
            if isempty(find(isinf(VARI{k}.data)))==0
                VARI{k}.data(isinf(VARI{k}.data))=NaN;
                warning(['infinite value found in ' VARI{k}.name ' :stored as missing_value in the netcdf file']);
            end
            % s'il n'y a  pas de missval, on en cree une
            % if isempty(VARI{k}.data(isnan(VARI{k}.data)))==0
            %if isempty(missval)
            %    missval=99999;
            %end
            % end

            dataTYP='float';
            if isa(VARI{k}.data,'double');dataTYP='double';end;
            if isinteger(VARI{k}.data);dataTYPe='int';end;
            if isfield(VARI{k},'type')
                dataTYP=VARI{k}.type;
            end

        end

        if ischar(VARI{k}.data)
            dataTYP='char';
            %if isempty(missval)
            %   missval=' ';
            %end
        end
    end

   
    % Definition  des dimensions
    % ------------------------------------------------------------------

    if k<=NBVAR1
        if isfield(VARI{k},'dimlength')==1
            expr{1}=['f(''' VARI{k}.name ''') = VARI{k}.dimlength;'];
            if strcmp(VARI{k}.name,REC)==1
            expr{1}=['f(''' VARI{k}.name ''') = 0;']; % record dimension
            end
        end
    end
    
    
    % ECRITURE  des variables
    % ------------------------------------------------------------------
    if isfield(VARI{k},'data')==1
        
        % definition des variables
        
        expr{2}=define_var(VARI{k},extract_dim,dataTYP);

        % recherche tous les attributs et leur type
        attnames=fieldnames(VARI{k});
        nbatt= length(attnames);
        jatt=2;
        for iatt=1:nbatt
            if verbose==1
                disp('**********************')
                disp(['attributs: ' attnames{iatt}])
                disp('-----------------')
            end
            theattname = attnames{iatt};

            if ~(strcmp(theattname,'data')==1|strcmp(theattname,'dim')==1|strcmp(theattname,'type')==1|strcmp(theattname,'name')==1)
                thevalatt = VARI{k}.(attnames{iatt});
                if verbose==1
                    thevalatt
                end
                
                if isempty(thevalatt)==0
                    jatt=jatt+1;
                    if ischar(thevalatt)
                        %disp('ok char')
                        expr{jatt}=['f{''' VARI{k}.name '''}.' theattname '=ncchar(''' thevalatt ''');'];
                        if verbose==1
                            jatt
                            expr{jatt}
                        end
                    elseif isfloat(thevalatt)
                        expr{jatt}=['f{''' VARI{k}.name '''}.' theattname '=ncfloat(' num2str(thevalatt) ');'];
                        if isa(thevalatt,'double')
                            expr{jatt}=['f{''' VARI{k}.name '''}.' theattname '=ncdouble(' num2str(thevalatt) ');'];
                        end
                        if verbose==1
                            jatt
                            expr{jatt}
                        end
                    elseif islogical(thevalatt)
                        expr{jatt}=['f{''' VARI{k}.name '''}.' theattname '=ncbyte(' num2str(thevalatt) ');'];
                    end
                    if strcmp(theattname,'valid_min')==1|strcmp(theattname,'valid_max')==1|strcmp(theattname,'resolution')==1
                        expr{jatt}=['f{''' VARI{k}.name '''}.' theattname '=nc' dataTYP '(' num2str(thevalatt) ');'];
                    end
                    expr{jatt};
                end
            end

        end
       
        if verbose==1
            pause
        end
        
        % store the data
        if isempty(VARI{k}.data(:))==0
            expr{jatt+1}=store_var(VARI{k},extract_dim,dataTYP,Dim,REC,REC_length);
            %disp(expr{jatt+1})
            %expr{jatt+1}=['f{''' VARI{k}.name '''}(:) = VARI{k}.data(:);'];
        end
       % pause
    end

    for ek=1:length(expr)
        if isempty(expr{ek})==0
          
            
            %expr{ek}
            %pause
            if verbose==1
                expr{ek}
            end
            %expr{ek}
            %VARI{k}.dimlength
            eval(expr{ek});
        end
    end

    if verbose==1
        disp('====================================== ')
    end

end


f=close(f);
% disp('fermeture du fichier')


%________________________________________________________________________________________________________________________
%________________________________________________________________________________________________________________________
%% SUBROUTINE
function [expression]=store_var(VAL,extract_dim,dataTYP,Dim,REC,REC_length)

% expr{jatt+1}=['f{''' VARI{k}.name '''}(:) = VARI{k}.data(:);'];
% si une des dimension et la "record dimension, alors il faut ecrire:
% f{''' VARI{k}.name '''}(1:dimrec,:,:) = VARI{k}.data(:) par exemple

siz=length(extract_dim);

fulldim=extract_dim;
thestr='(';
for ndim=1:siz
    
    %thestr = [thestr '1:' num2str(Dim.(lower(fulldim{ndim})).dimlength) ','];
     thestr = [thestr '1:' num2str(Dim.(lower(fulldim{ndim})).dimlength) ','];
end
thestr=[thestr,')'];
thestr=strrep(thestr,',)',')');


expression=['f{''' VAL.name '''}' thestr '= VARI{k}.data(:);'];


return
%________________________________________________________________________________________________________________________

function [expression]=define_var(VAL,extract_dim,dataTYP)

siz=length(extract_dim);

fulldim=extract_dim;

if(siz>=1)
    expression=['f{''' VAL.name '''} = nc' dataTYP '(''' fulldim{1} ''''];
end

for k=2:siz
    expression=[expression ',''' fulldim{k} ''''];
end

expression=[expression ');'];


return

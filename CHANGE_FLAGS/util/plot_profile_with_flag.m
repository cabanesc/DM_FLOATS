% -========================================================
%   USAGE : [thetitle]=plot_profile_with_flag(S,ParamX,ParamY,profile_number,varargin)
%   PURPOSE : plot d'un profil PSAL ou TEMP avec les flags associes
% -----------------------------------
%   INPUT :
%     S   (structure)  - float data (multiprofile or single)
%     ParamX(char)     - parametre à tracer ('psal' ou 'temp' ou 'tpot' ou 'pres')
%     ParamY(char)     - parametre à tracer ('psal' ou 'temp' ou 'tpot' ou 'pres') )
%     profile_number  (integer) numero de profil (peut etre different de numero de cycle
%   OPTIONAL INPUT
%     '
%     'LineWidth' (integer) epaisseur des lignes
%     'LineColor' (character) couleur des lignes (parmis 'b','r','g','m','c','k','w')
%     'Marker'    (character) matlab markers ('+','.','^' or 'none')
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx
%   CALLED SUBROUTINES: none
% ========================================================

function [thetitle]=plot_profile_with_flag(S,ParamX,ParamY,profile_number,varargin)



n=length(varargin);
p=struct;

% verifie que la premiere dimension de tous les tableaux est N_PROF (rearrange sinon)
S = check_FirstDimArray_is(S,'N_PROF');
map=jet(S.n_prof);
if n>0
    if n/2~=floor(n/2)
        error('check the imput arguments')
    end
    
    f=varargin(1:2:end);
    c=varargin(2:2:end);
    p = cell2struct(c,f,2);
end

if isfield(p,'Marker')==0
    p.Marker='+';
end
if isfield(p,'LineWidth')==0
    p.LineWidth=3;
end
if isfield(S,'max_cy')&&isfield(p,'LineColor')==0
Vec=[0:S.max_cy];
c=colormap(jet(length(Vec)));
p.LineColor='multi';
end
if isfield(p,'LineColor')==0
    p.LineColor='b';
end

if isfield (S, 'fillisnan')
    if S.fillisnan==0
        S = replace_fill_bynan(S);
    end
else
    S = replace_fill_bynan(S);
end
thetitle=[];

ParamX = lower(ParamX);
ParamY = lower(ParamY);

nn=length(profile_number);
for ikk=1:nn
    nprof=profile_number(ikk);
    
    if strcmp(p.LineColor,'multi')
       ii=find(Vec==S.cycle_number.data(nprof));
       theLineColor=c(ii,:);
    else
    theLineColor=p.LineColor;
    end
    if isfield(S,ParamX)  % abscisse
        if isfield(S.(ParamX),'data')
            if nprof>0&nprof<=size(S.(ParamX).data,1)
                subX=S.(ParamX).data(nprof,:);
                if isfield(S,ParamY)  % abscisse
                    if isfield(S.(ParamY),'data') 
						
						level = S.(ParamY).data(nprof,:);
                        plot(subX,level ,'Color',theLineColor,'LineWidth',p.LineWidth) ;
                        
                        hold on
                        grid on
						
                        xlabel([upper(ParamX)],'interpreter','none')
                        ylabel([upper(ParamY)],'interpreter','none')
                        set(gca,'Ydir','normal')
                        if isempty(findstr(ParamY,'pres'))==0
                            set(gca,'Ydir','reverse')
                        end
                        
                        if isfield(S, [ParamX '_qc'])
                            if isfield(S.([ParamX '_qc']),'ischar2num')
                                if S.([ParamX '_qc']).ischar2num~=1
                                    warning('Flag values are characters, should be numeric for plotting purpose, you can use format_flags_char2num.m ')
                                end
                            else
                                %  warning('Check if flag values are numeric, you can use format_flags_char2num.m ')
                            end
                            
                            if isfield(S.([ParamX '_qc']),'data')
							
                                the_qc = S.([ParamX '_qc']).data(nprof,:);   
                                sub2=subX;
                                sub2(the_qc~=1)=NaN;
								%p.Marker
                                if strcmp(p.Marker,'none')==0
                                    
                                    colorp=['g' p.Marker];
                                   
                                    plot(sub2,level ,colorp,'LineWidth',p.LineWidth) ;
                                    
                                    colorp=['y' p.Marker];
                                    sub2=subX;
                                    sub2(the_qc~=2)=NaN;
                                    plot(sub2,level ,colorp,'LineWidth',p.LineWidth) ;
                                    
                                    colorp=['m' p.Marker];
                                    sub2=subX;
                                    sub2(the_qc~=3)=NaN;
                                    plot(sub2,level ,colorp,'LineWidth',p.LineWidth*2) ;
                                    
                                    colorp=['r' p.Marker];
                                    sub2=subX;
                                    sub2(the_qc~=4)=NaN;
                                    plot(sub2,level ,colorp,'LineWidth',p.LineWidth*2) ;
                                    
                                    colorp=['b' p.Marker];
                                    sub2=subX;
                                    sub2(the_qc==1|the_qc==2|the_qc==3|the_qc==4)=NaN;
                                    plot(sub2,level ,colorp,'LineWidth',p.LineWidth) ;
                                end
                            end
                        end
                        
                    end
                    
                end
            end
        end
    end
    
    
    % titre du plot (numero de flotteur, date, longitude, latitude, cycle, ) si on a les infos
    
    thetitle = [' '];
    info1  = {'platform_number','juld','longitude','latitude','cycle_number','direction','pi_name','data_centre'};
    info2  = {'', ', ', ', lon: ',', lat: ',', cycle: ' ,'', ', PI: ', ', '};
    
%    S.cycle_number.data(nprof,:)=num2str(S.cycle_number.data(nprof,:));
    
    if nprof>0&nprof<=size(S.(ParamX).data,1)
        %keyboard
        for k=1:length(info1)
            if isfield(S,info1{k})
                if isfield(S.(info1{k}), 'data')
                    if isequal(info1{k},'juld')
                        % transforme la date juld > jour dd/mm/yyyy
                        thedate = datestr( (S.juld.data(nprof,:) + datenum('19500101','yyyymmdd')),'dd/mm/yyyy');
                        
                        thetitle = [thetitle info2{k}  thedate  ];
                    elseif isequal(info1{k},'cycle_number')
                        if nn>1
                            thecycles=[num2str(S.cycle_number.data(1)) '-' num2str(S.cycle_number.data(end))];
                            info2{k}=', cycles: ';
                        else
                            thecycles=[num2str(S.cycle_number.data(nprof))];
                        end
                        thetitle = [thetitle info2{k}  thecycles  ];
                    else
                        thetitle = [thetitle info2{k}  strtrim(num2str(S.(info1{k}).data(nprof,:)))  ];
                    end
                    
                end
            end
        end
    else
        thetitle='';
    end
    
    %title(thetitle,'interpreter','none')
end
thetitle
return

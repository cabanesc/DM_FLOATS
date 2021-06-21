% -========================================================
%   USAGE : [thetitle]=plot_profile(S,Param,profile_number,varargin)
%   PURPOSE : plot d'un profil PSAL ou TEMP 
% -----------------------------------
%   INPUT :
%     S   (structure)  - float data (multiprofile or single)
%     Param (char)     - parametre à tracer ('psal' ou 'temp')
%     profile_number  (integer) numero de profil (peut etre different de numero de cycle
%   OPTIONAL INPUT
%     'Mindepth'  (double) profondeur minimum du plot
%     'LineWidth' (integer) epaisseur des lignes
%     'LineColor' (character) couleur des lignes (parmis 'b','r','g','m','c','k','w')
%     'Marker'    (character) matlab markers ('+','.','^' or 'none')
%     'OnlyGood'  (character) 'y' or 'n'. If yes keep only data with flags <=2
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx
%   CALLED SUBROUTINES: none
% ========================================================

function [thetitle]=plot_profile(S,ParamX,ParamY,profile_number,varargin)


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
if isfield(p,'LineColor')==0
   %p.LineColor='b';
end
if isfield(p,'Marker')==0
   p.Marker='none';
end
if isfield(p,'LineWidth')==0
   p.LineWidth=2;
end
if isfield(p,'OnlyGood')==0
   p.OnlyGood='y';
end

S = format_flags_char2num(S);
if strcmp(p.OnlyGood,'y')
S.psal.data(S.psal_qc.data>2)=NaN;
S.temp.data(S.temp_qc.data>2)=NaN;
end

if ~isfield(S,'tpot')
       S.tpot.data = sw_ptmp(S.psal.data,S.temp.data,S.pres.data,1000);
else
       S.tpot.data(S.temp_qc.data>2|S.psal_qc.data>2)=NaN;
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
nprof=profile_number;
%keyboard

if isfield(S,ParamX)  % abscisse
    if isfield(S.(ParamX),'data')
        if profile_number>0&profile_number<=size(S.(ParamX).data,1)
            subX=S.(ParamX).data(nprof,:);
        
        
            if isfield(S,ParamY)  % abscisse
                if isfield(S.(ParamY),'data')
            
                
                    level = S.(ParamY).data(nprof,:);
                    
                    %keyboard
                    
                    if isfield(p,'LineColor')==1
                    for ik=1:length(nprof)
                        if isfield(p,'Mindepth')==1
                            pm=find(level(ik,:)>p.Mindepth);
                            plot(subX(ik,pm),level(ik,pm) ,'Color',p.LineColor,'LineWidth',p.LineWidth,'Marker',p.Marker) ;
                        else
                            plot(subX(ik,:),level(ik,:) ,'Color',p.LineColor,'LineWidth',p.LineWidth,'Marker',p.Marker) ;
                        end
                    end
                    else
                        
                      for ik=1:length(nprof)
                        if isfield(p,'Mindepth')==1
                            pm=find(level(ik,:)>p.Mindepth);  
                            plot(subX(ik,pm),level(ik,pm) ,'Color',map(nprof(ik),:),'LineWidth',p.LineWidth,'Marker',p.Marker) ;
                        else
                        plot(subX(ik,:),level(ik,:) ,'Color',map(nprof(ik),:),'LineWidth',p.LineWidth,'Marker',p.Marker) ;
                        end
                      end
                    end
                    hold on
                    grid on
                    
                    xlabel([upper(ParamX)],'interpreter','none')
                    ylabel([upper(ParamY)],'interpreter','none')
                    
                    if isempty(findstr(ParamY,'pres'))==0
                        set(gca,'Ydir','reverse')
                    end
                    

                end
                
                
            end
        end
    end
    
    % titre du plot (numero de flotteur, date, longitude, latitude, cycle, ) si on a les infos
    thetitle = [strtrim(num2str(S.platform_number.data(1,:))),', cycles:' strjoin(cellstr(num2str(S.cycle_number.data(profile_number)'))',',')] ;
    
    
end
box on
grid on

return

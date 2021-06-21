% -========================================================
%   USAGE : [thetitle]=plot_theta_S_diag(S,profile_number,isopycne,varargin)
%   PURPOSE : plot du diagramme tetha-S d'un ou plusieurs profils
% -----------------------------------
%   INPUT :
%     S   (structure)  - float data (multiprofile or single)
%     profile_number  (integer) numero de profil (peut etre different de numero de cycle)
%     isopycne (character) 'y' or 'n' pour tracer ou pas les isopycnes sur le diagramme
%   OPTIONAL INPUT
%     'LineWidth' (integer) epaisseur des lignes
%     'LineColor' (character) couleur des lignes (parmis 'b','r','g','m','c','k','w')
%     'Marker'    (character) matlab markers ('+','.','^' or 'none')
%     'OnlyGood'  (character) 'y' or 'n'. If yes keep only data with flags <=2
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx
%   CALLED SUBROUTINES: none
% ========================================================

function [thetitle]=plot_theta_S_diag(S,profile_number,isopycne,varargin)


n=length(varargin);
p=struct;

% verifie que la premiere dimension de tous les tableaux est N_PROF (rearrange sinon)
S = check_FirstDimArray_is(S,'N_PROF');
if size(profile_number,1)>1
profile_number=profile_number';
end

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
  % p.LineColor='b';
end
if isfield(p,'Marker')==0
   p.Marker='.';
end
if isfield(p,'LineWidth')==0
   p.LineWidth=2;
end
if isfield(p,'OnlyGood')==0
   p.OnlyGood='y';
end


if isfield (S, 'fillisnan')
    if S.fillisnan==0
        S = replace_fill_bynan(S);
    end
else
    S = replace_fill_bynan(S);
end
thetitle=[];

S = format_flags_char2num(S);
if strcmp(p.OnlyGood,'y')
S.psal.data(S.psal_qc.data>2)=NaN;
S.temp.data(S.temp_qc.data>2)=NaN;
end

%

if isfield(S,'psal')& isfield(S,'temp')&isfield(S,'pres')
   if ~isfield(S,'tpot')
       S.tpot.data = sw_ptmp(S.psal.data,S.temp.data,S.pres.data,1000);
   else
      % S.tpot.data(S.temp_qc.data>2|S.psal_qc.data>2)=NaN;
   end
   for iprof=1:length(profile_number)
      if isfield(p,'LineColor')==0
      plot(S.psal.data(profile_number(iprof),:),S.tpot.data(profile_number(iprof),:),'color',map(profile_number(iprof),:),'marker',p.Marker,'linewidth',p.LineWidth);
      else
      plot(S.psal.data(profile_number(iprof),:),S.tpot.data(profile_number(iprof),:),'color',p.LineColor,'marker',p.Marker,'linewidth',p.LineWidth);
      end
    end
    if isopycne=='y'
       ctpot=[floor(min(min(S.tpot.data(profile_number,:))))-1:0.1:ceil(max(max((S.tpot.data(profile_number,:)))))+1]';
       cpsal=[floor(10*min(min(S.psal.data(profile_number,:))))/10-0.1:0.1:ceil(max(max(S.psal.data(profile_number,:)))*10)/10+0.1]';
       [tabct,tabcp]=meshgrid(ctpot,cpsal);
       [~,csig]=swstat90(tabcp,tabct,0);
       if size(csig,1)~=size(tabcp,1)
          csig=csig';
       end
       if sum(isnan(tabcp))~=length(tabcp)& sum(isnan(tabct))~=length(tabct)
       [c,h]=contour(tabcp,tabct,csig,[20:0.5:35],':k');
       clabel(c,h);
       end
    end
    ylabel('Potential Temp. (ref. to 0db)')
    xlabel('Salinity')
    vxmax=max(max(ceil(S.psal.data(profile_number,:)*100)))/100;
    vxmin=min(min(floor(S.psal.data(profile_number,:)*100)))/100;
    vymax=max(max(ceil(S.tpot.data(profile_number,:)*100)))/100;
    vymin=min(min(floor(S.tpot.data(profile_number,:)*100)))/100;
    if ~isnan(vxmin+vxmax+vymin+vymax)
    set(gca,'xlim',[vxmin vxmax],'ylim',[vymin vymax]);
    end
    grid on
    box on 

    
    % titre du plot (numero de flotteur, date, longitude, latitude, cycle, ) si on a les infos
    
    thetitle = [strtrim(num2str(S.platform_number.data(1,:))),', cycles:' strjoin(cellstr(num2str(profile_number'))',',')] ;
    
    
    %title(thetitle,'interpreter','none')
    
end

return

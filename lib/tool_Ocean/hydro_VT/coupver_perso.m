function [M2,x2,p2,a1,a2,Cca,hca,Ccaf,hcaf,tdate,tref]=coupver(fichent,param,sta,varargin)
%function [M,x,p,a1,a2,Cca,hca,tdate,tref]=coupver(fichent,param,sta)
%M=matrice totale des valeurs du parametre choisi
%x,p : absisses et pressions de l'ensemble de la section
%a1 : vecteur des axes (sens p croissant) de la derniere figure affichee 
%a2 : vecteur des axes globaux (titre et textes)
%fichent=nom du fichier d'entree netcdf multi-station dans /home/menbrial/HYDRO_ATLN/MLT_NC/ 
%        (de 1 en 1 metre)
%param=nom du parametre a tracer (identique a celui dans le fichier). Ex.:
%      'TEMP'  temperature
%      'TPOT' temperature potentielle
%      'PSAL'  salinite
%      'SIG0'  densite g(0,tp,s)
%      'SIGI' densite g(p,t,s)
%      'OXYL' oxygene (ml/l)
%      'DYNH' hauteur dynamique (ref. a 0m par defaut)     %%%
%      'VGEO'  vitesse geostrophique perpendiculaire a la section  %%%
%      'BRV2'  Periode de Brunt-Vaisala (heures)  %N%
%      'VORP'  Vorticite Potentielle (m-1s-1) %-(N2*f/g)
%sta=vecteur des stations de la section
%
%En optionnel :
%'color'			: section en couleurs degradees
%'titre','...' 			: titre du dessin (en haut)
%'comment','...'		: commentaire (en bas)
%'portrait'                     : format portrait (paysage par def.)
%'position',[left bot wid h]    : position de la fig (cm) dans la page 
%                                 pleine page par defaut (nan=>defaut)
%'hold'                         : ajoute la section sur la figure existante
%'degres','dec' ou 'min'	: degres decimaux par defaut
%'plim',[pmin p1 p2 ... pmax]	: les limites en pression
%'dpression',dp			: ecart entre les pressions considerees
%				  pour les contours (10dbar par defaut)
%'pourcentage',[po1 ...pon]	: proportion de chaque partie en %
%'proj',proj			:    proj.='long' (par defaut)
%				  ou proj.='lat' ou proj='mille' ou 'km'
%'xlim',[xmin xmax]		: limites abscisses (unites de la projection)
%'ech',ech			: nombre de cm par deg ou par 100 milles ou km;
%				  max=25cm au total; par defaut : pleine ech
%'dpa',[dpa1 dpa2 ...]		: ecart entre les isolignes dans chaque partie
%'vpa',[pa11 pa12 ...;pa21 ...]	: vecteur des valeurs des isolignes par partie
%'vpb',[pa11,pa12 ...;pa21 ...]	: vecteur des valeurs des isolignes en gras
%'contour','off'		: n'affiche pas les contours ('on' par defaut)
%'clab','off'       : n'affiche pas les labels des contours ('on' par defaut)
%'matrice','on'			: affiche la matrice M ('off' par defaut)
%'pasmat',pasmat		: pas de l'ecriture sur y (10db par defaut)
%				  et passe 'matrice' a la valeur 1
%'ref',ref			: pression reference pour hd ou u (fond (inf) par defaut).	%%%
%'opt',iopt                     : choisir 1 des 4 solutions pour les bords 
%                                 iopt=1 par defaut. iopt=0 => 4 solutions affichees. 
%'fill',[v11 v12 c1; v21 v22...]: colore en c1 entre les contours v11 et v12,
%				  en c2 entre les contours v21 et v22 ...
%				  v peut avoir comme valeur inf et -inf
%				  c1,c2... sont des vecteurs lignes RGB
%'fill',[v1 v2]                 : colore entre tous les contours fins (eviter les
%                                 gras dans cette option) en ajustant la palette
%                                 clr1 entre les valeurs v1 et v2. NaN pour v1 ou
%                                 v2 => ajustement sur les val. ext. calculees.
%'liss',lis     : nombre de metres de lissage (50m par defaut)

% Modifs : aller-retour sur x pour completer les 'trous'
%	   jeg pour le calcul de gradvia
%	   vpa,vpb (dpa):1 seule ligne (valeur) acceptee meme si plusieurs axes
%	   xlim introduit

%exemples pour OVIDE 2002:
%
%load mapcolor2
%coupver('OVIDE/ovid02_dep.nc','TPOT',[6:72 74:95],'vpa',[0:0.25:3,3.5:0.5:7,8:17],'fill',[1 16],'tit','OVIDE 2002','com','Potential Temperature','proj','km');
%colormap(jetemp);
%coupver('OVIDE/ovid02_dep.nc','PSAL',[6:72 74:95],'vpa',[33 34.800:0.02:35,35.1:0.1:36.3],'fill',[34.82 36.2],'tit','OVIDE 2002','com','Salinity','proj','km');
%colormap(jetsal);
%coupver('OVIDE/ovid02_dep.nc','OXYL',[6:72 74:95],'vpa',[4:.25:9],'fill',[4 7],'tit','OVIDE 2002','com','Oxygen (ml/l)','proj','km');
%colormap(vmap)
%coupver('OVIDE/ovid02_dep.nc','VGEO',[6:72 74:95],'fill',[-40 40],'vpa',[-40:10:60],'tit','OVIDE 2002','com','Geostrophic velovity (cm/s)','proj','km');
%colormap(vmap)
%coupver('OVIDE/ovid02_dep.nc','SIG0',[6:72 74:95],'fill',[26 27.92],'vpa',[26:.5:27 27.2:.2:27.6 27.7:.05:27.85 27.88:.01:27.9],'tit','OVIDE 2002','com','Density anomaly referenced to 0m','proj','km');
%coupver('OVIDE/ovid02_dep.nc','SIG1',[6:72 74:95],'fill',[30.2 32.6],'vpa',[30 31 31.2:.2:31.8 32:.1:32.4 32.45:.05:32.65],'tit','OVIDE 2002','com','Density anomaly referenced to 1000m','proj','km');
%coupver('OVIDE/ovid02_dep.nc','SIG2',[6:72 74:95],'fill',[34.5 37.2],'vpa',[34.5:.5:36.5 36.6:.1:36.9 36.95:.05:37.2],'tit','OVIDE 2002','com','Density anomaly referenced to 2000m','proj','km');
%coupver('OVIDE/ovid02_dep.nc','BRV2',[6:72 74:95],'fill',[0 5],'vpa',[0 0.5 1:4 8 40],'tit','OVIDE 2002','com','Buoyancy period (hour)','proj','km');
%colormap(jetpers);
%coupver('OVIDE/ovid02_dep.nc','N2',[6:72 74:95],'fill',[0 5],'vpa',[0 0.5 1:4 8 40],'tit','OVIDE 2002','com','Buoyancy period (hour)','proj','km');
%colormap(jetpers);
%coupver('OVIDE/ovid02_dep.nc','VORP',[6:72 74:95],'fill',[-0.01 0.2],'vpa',[0:0.01:0.02 0.05:0.05:0.2 1],'tit','OVIDE 2002','com','Potential vorticity (10^{-9} s^{-1})','proj','km');
%colormap(jetpers);

% Filtrage sur lis valeurs espacees de dp.

%%initialisation (pour eviter les warning dus aux isempty)
M=[];
titre=[]; com=[]; dpa=[]; vpa=[]; vpb=[]; plim=[]; dp=10; pourc=[]; proj='lo';
xlim=[]; ech=[]; cont='on'; contlab='on'; mat=[]; dmat=[];deg='dec'; iopt=1; ref=inf; fil=[];
col=0; ori='pay'; posi=[]; holdon=logical(0); lis=50;

format compact
%% traitememt des entrees et etablissement des valeurs par defaut
ii=0;
while (ii<=nargin-4),
  ii=ii+1;
  arg=varargin{ii};
  if (isstr(arg)),
    if     lower(arg(1:3))=='col', col=1;
    elseif lower(arg(1:3))=='tit', ii=ii+1; titre=varargin{ii};
    elseif lower(arg(1:3))=='com', ii=ii+1; com=varargin{ii};
    elseif lower(arg(1:3))=='dpa', ii=ii+1; dpa=varargin{ii};
    elseif lower(arg(1:3))=='vpa', ii=ii+1; vpa=(varargin{ii})';
    elseif lower(arg(1:3))=='vpb', ii=ii+1; vpb=(varargin{ii})';
    elseif lower(arg(1:3))=='pli', ii=ii+1; plim=varargin{ii};
    elseif lower(arg(1:3))=='dpr', ii=ii+1; dp=varargin{ii};
    elseif lower(arg(1:3))=='pou', ii=ii+1; pourc=varargin{ii};
    elseif lower(arg(1:3))=='lis', ii=ii+1; lis=varargin{ii};    
    elseif lower(arg(1:3))=='pro', ii=ii+1; proj=varargin{ii};
    elseif lower(arg(1:3))=='xli', ii=ii+1; xlim=varargin{ii};
    elseif lower(arg(1:3))=='ech', ii=ii+1; ech=varargin{ii};
    elseif lower(arg(1:3))=='con', ii=ii+1; cont=varargin{ii};
    elseif lower(arg(1:3))=='cla', ii=ii+1; contlab=varargin{ii};
    elseif lower(arg(1:3))=='mat', ii=ii+1; mat=varargin{ii};
    elseif lower(arg(1:3))=='pas', ii=ii+1; dmat=varargin{ii};
    elseif lower(arg(1:3))=='deg', ii=ii+1; deg=varargin{ii};
    elseif lower(arg(1:3))=='opt', ii=ii+1; iopt=varargin{ii};
    elseif lower(arg(1:3))=='ref', ii=ii+1; ref=varargin{ii};  %%%
    elseif lower(arg(1:3))=='fil', ii=ii+1; fil=varargin{ii};    
    elseif lower(arg(1:3))=='por', ori='por';    
    elseif lower(arg(1:3))=='pay', ori='pay';    
    elseif lower(arg(1:3))=='pos', ii=ii+1; posi=varargin{ii};    
    elseif lower(arg(1:3))=='hol', holdon=logical(1);    
    else error([arg ' : Parametre inconnu']);
    end;
  end;
end;
para=param;
if strcmp(param,'VGEO'),  para='DYNH'; end;
if strcmp(param,'N2'),  para='TEMP'; end;
%if strcmp(param,'VPOT'),  para='BRV2'; end;

%if     strcmp(param,'u'),  npa=3;
%elseif strcmp(param,'tp'), npa=4;
%elseif strcmp(param,'s'),  npa=5;
%elseif strcmp(param,'g'),  npa=6;  %en fait, calcul refait ici
%elseif strcmp(param,'gt'), npa=7;
%elseif strcmp(param,'h'), npa=8;   %en fait, calcul refait ici
%elseif strcmp(param,'o2'), npa=10;
%elseif strcmp(param,'t'), npa=3;   %calcul fait ici
%elseif strcmp(param,'N'), npa=6;   %calcul fait ici avec t et s (Nbv2.m)%N%
%end;
nsta=length(sta);
if nsta<2, error('Au moins 2 stations !'); end;
proj=lower(proj(1:2));
if sum(proj=='mi' | proj=='lo' | proj=='la' | proj=='km')'~=2,
  error('La projection n''est pas correctement definie');
end;
if ~isempty(ech) & (proj=='mi' | proj=='km'), ech=ech/100; end;
if isempty(mat) & isempty(dmat), mat='off';
elseif isempty(dmat), dmat=10; 
elseif isempty(mat), mat='on';
end;
if ~isempty(plim),
  plim=sort(plim);
  lplim=length(plim);
  if lplim>2
    if ~isempty(vpa), lvpa=size(vpa,2);
                      if lvpa==1, vpa=vpa*ones(1,lplim-1);
                      elseif lvpa~=lplim-1, error('vpa incorrect'); 
                      end;
    end;
    if ~isempty(vpb), if size(vpb,2)==1, vpb=vpb*ones(1,lplim-1);
                      elseif size(vpb,2)~=lplim-1, error('vpb incorrect');
                      end;
    end;
    if ~isempty(dpa), if length(dpa)==1, dpa=ones(lplim-1,1)*dpa;
                      elseif length(dpa)~=lplim-1, error('dpa incorrect');
                      end;
    end;
    if isempty(pourc), pourc=100*diff(plim)/(max(plim)-min(plim)); end;
  elseif lplim==2, pourc=100;
  elseif lplim<2, 
         error('La propriete plim n''est pas correctement definie');
  end;
end;
if ~isempty(pourc) 
 if isempty(plim), error('Vous devez definir des pressions limites');
 elseif length(pourc)~=length(plim)-1, error('erreur le nombre de pourcentages');
 end;
end;
lef=500; bot=50; larg=750; haut=950;
if iopt==0, iopt=1:4; lef=250; bot=350; larg=250; haut=600; end;

if strcmp(ori,'pay'), Lfig=27.7; Hfig=19; else Lfig=19; Hfig=27.7; end;
if isempty(posi),
  posi=[.1*Lfig .15*Hfig .8*Lfig .6*Hfig];
  if ~isempty(ech), posi(1)=nan; end;
else 
  if isnan(posi(3)),
    if isnan(posi(1)), posi(3)=.8*Lfig; else posi(3)=.9*Lfig-posi(1); end;
  end;
  if isnan(posi(2)), posi(2)=.13*Hfig; end;
  if isnan(posi(4)), posi(4)=.845*Hfig-posi(2); end;
end;
if posi(2)+posi(4)>.9*Hfig,
 error(sprintf('La figure ne fait que %4.1f cm de hauteur\n',Hfig*.9)); 
end;
if posi(1)+posi(3)>.9*Lfig,
 error(sprintf('La figure ne fait que %4.1 cm de largeur\n',Lfig*.9)); 
end;

%% chargement des donnees
if strcmp(fichent(1),'/'),
  direc='';
else
  direc='/home1/penfret/HYDROCEAN/MLT_NC/';
end
if strcmp(fichent(end-1:end),'nc'),
  ncload([direc fichent],para,'DEPH','PRES','BOTTOM_DEPTH','MAX_VALUE_PARAM_REF',...
         'LATITUDE_BEGIN','LATITUDE_END','LONGITUDE_BEGIN','LONGITUDE_END',...
         'STATION_NUMBER','JULD_BEGIN','JULD_END','MAX_PRESSURE');
elseif strcmp(fichent(end-2:end),'mat'),
  load([direc fichent]);
else
  error('Unknown file extension');
end
     
DEPH(DEPH==-9999)=NaN;
PRES(PRES==-9999)=NaN;
para=eval(para);
para(para==-9999)=NaN;
LATITUDE_BEGIN(LATITUDE_BEGIN==-9999)=NaN;
LATITUDE_END(LATITUDE_END==-9999)=NaN;
LONGITUDE_BEGIN(LONGITUDE_BEGIN==-9999)=NaN;
LONGITUDE_END(LONGITUDE_END==-9999)=NaN;
JULD_BEGIN(JULD_BEGIN==-9999)=NaN;
JULD_END(JULD_END==-9999)=NaN;
LATITUDE=mynanmean([LATITUDE_BEGIN';LATITUDE_END'])';
LONGITUDE=mynanmean([LONGITUDE_BEGIN';LONGITUDE_END'])';
JULD=mynanmean([JULD_BEGIN';JULD_END'])';
%%extremes de p (on garde le nom p, mais on raisonne en profondeur)
%for i=1:nsta,
% ista=find(STATION_NUMBER==sta(i));
% %a=fichent(fichent(:,1)==sta(i),:);
% if isempty(ista), 
%  fprintf('La station %g n''existe pas dans le fichier\n',sta(i));
% else
%  %p=a(2:size(a,1),2);
%  if i==1, pmin=min(p); pmax=max(p);
%  else if min(p)<pmin, pmin=min(p); end;
%       if max(p)>pmax, pmax=max(p); end;
%  end;
% end;
%end;
pmin=min(DEPH(:));
pmax=max(DEPH(:));
pmax=pmax+2*dp;
if isempty(plim), plim=[pmin pmax]; ledge=2; pourc=100; end;
%p=(pmin:dp:pmax)';
p=(plim(1):dp:plim(end))';
lp=length(p);
%
%%boucle de construction de la matrice des parametres M
ista=0; %ista indice les stations qui satisfont aux conditions de l'utilisateur
for i=1:nsta,
  ii=find(STATION_NUMBER==sta(i));
  %a=fichent(fichent(:,1)==sta(i),:);
  if isempty(ii), 
    b=[]; 
    fprintf('La station %g n''existe pas dans le fichier\n',sta(i));
  else 
    isok=find(~isnan(para(ii,:)));
    b=[DEPH(ii,isok)' para(ii,isok)'];
    iplim=find(b(:,1)>=pmin & b(:,1)<=pmax);
    b=b(iplim,:);
  end;
  if ~isempty(b),
   ista=ista+1;
   staf(ista)=sta(i); 
   lat(ista)=LATITUDE(ii); 
   lg(ista)=LONGITUDE(ii);
%   lso(ista)=profpres(a(1,8),lat(ista)); %profondeur de la station en dbar
   lso(ista)=BOTTOM_DEPTH(ii); %profondeur de la station en m
   max_pres(ista)=MAX_PRESSURE(ii); %max pression des mesures
   refhd(ista)=ref(1);
   dat=JULD(ii);
   if ista==1, jjmi=dat; jjma=dat;
   elseif dat<jjmi, jjmi=dat;
   elseif dat>jjma, jjma=dat;
   end
%   dat=a(1,2:4);
%   if ista==1, ami=dat(3); ama=dat(3); jjmi=j2jul(dat); jjma=jjmi;
%   else
%      if dat(3)>ama, ama=dat(3); jjma=j2jul(dat); end;
%      if dat(3)<ami, ami=dat(3); jjmi=j2jul(dat); end;
%      if (dat(3)==ama & j2jul(dat)>jjma), jjma=j2jul(dat); end;
%      if (dat(3)==ami & j2jul(dat)<jjmi), jjmi=j2jul(dat); end;
%   end;
   ps=b(:,1); par=b(:,2);
   if strcmp(param,'BRV2'),
      par(par<2.1154e-8)=2.1154e-8;  %%eq a 12h
      par=2*pi./sqrt(par)/3600; 
      %par=moygliss(par,21,0,0);
   end;     %N%
   %if strcmp(param,'VPOT'),
      
   %end;     %N%
   if strcmp(param,'N2'), 
     if strcmp(fichent(end-1:end),'nc'),
       ncload([direc fichent],'PSAL');
     end
     b2=[PSAL(ii,isok)'];
     b2(b2==-9999)=NaN;
     b2=b2(iplim);
     [Nbv,par]=Nbv2(b2,b(:,2),b(:,1),21); 
   end;     %N%
   if strcmp(param,'VORP'),
      par=par*1e9;
   end
   if strcmp(param,'DYNH') | strcmp(param,'VGEO'),   %%%
      par=par*100;	%cmdyn
      if     length(ref)>1, ip=max(find(PRES(ii,isok)<=ref(i)));
      else   ip=max(find(PRES(ii,isok)<=ref));
      end
      refhd(ista)=PRES(ii,isok(ip));
      par=par-par(ip);
   end; %%%
   psmi=min(ps); psma=max(ps); 	%%bornes des valeurs experimentales
   par=spline(ps,par,p);	%% reechantillonnage selon dpr et plim
   par(p<psmi)=NaN*ones(size(par(p<psmi)));
   ipsm(ista)=max(find(p<=psma)); %% indice derniere valeur mesuree colonne i
   psm(ista)=p(ipsm(ista));	%% profondeur de la derniere valeur mesuree
   ifd(ista)=max(find(p<=lso(ista)));	%% indice du fond colonne i
%%  grad vertical moyen au fond :
   ia=min(ipsm(ista)-1,4);
   dparv(ista)=(par(ipsm(ista))-par(ipsm(ista)-ia))/ia; 
%% arrangement des valeurs au fond en vue de la moyenne
   if ipsm(ista)<lp-1,
    par(ipsm(ista)+1)=par(ipsm(ista))+dparv(ista);
    par(ipsm(ista)+2)=par(ipsm(ista)+1)+dparv(ista);%% + 2 valeurs au fond
    par(ipsm(ista)+3:lp)=NaN*ones(size(par(ipsm(ista)+3:lp)));
   elseif ipsm(ista)==lp-1, par(lp)=par(lp-1)+dparv(ista);
                            parfin=[par(lp-1);NaN];
   else parfin=[par(lp-1);par(lp)];
   end;
   %pardebut=par(1:2);
%   Mmat=[Mmat,par]; %% matrice avant filtrage pour affichage
%   par=moygliss(par,5,1); 
   %par(1:2)=pardebut; 
   %if ipsm(ista)>=lp-1, par(lp-1:lp)=parfin; end;
   M=[M,par]; %% matrice apres filtrage pour contours
  end;
end;
clear fichent
%filtrage final
if lis>=3*dp,
  M=filt_param(M',round(lis/dp))';
  fprintf('Les donnees sont decimees tous les %i dbars, puis filtrees sur %im\n',[dp lis]);
else
  fprintf('Les donnees sont decimees tous les %i dbars, mais pas filtrees\n',dp);
end

if ista<2, error('Il n''y a pas au moins 2 stations dans ces limites'); end;
%datmi=jul2j(jjmi,ami); datma=jul2j(jjma,ama);
datmi=greg_0h(jjmi+jul_0h(1950,1,1));
datma=greg_0h(jjma+jul_0h(1950,1,1));
dattext=sprintf('%i.%i,%i to %i.%i,%i',[datmi(2) datmi(3) datmi(1) ...
					datma(2) datma(3) datma(1)]);

%% Evaluation des abscisses x et reorganisation par ordre croissant de x
if strcmp(deg,'min'), lat=min2dec(lat); lg=min2dec(lg); end;
if proj=='mi'  | proj=='km',
  %islg=((max(lg)-min(lg))*cos(lat(1)/180*pi)>=(max(lat)-min(lat)));
  %if islg, [lg,ordre]=sort(lg'); lat=(lat(ordre))';
  %else [lat,ordre]=sort(lat'); lg=(lg(ordre))';
  %end; 
  ordre=[1:length(lg)]';
 
  dx=[0;(dist(lat,lg))'];
  x=cumsum(dx);
  if proj=='mi', x=x/1852; xlabe='distance (milles)';
  else x=x/1000; xlabe='Distance (km)';
  end;
  islg=2;
elseif proj=='lo', [lg,ordre]=sort(lg'); lat=(lat(ordre))'; x=lg; 
       islg=1; xlabe='longitude';
elseif proj=='la', [lat,ordre]=sort(lat'); lg=(lg(ordre))'; x=lat; 
       islg=0; xlabe='latitude';
end;
lso=(lso(ordre))';
staf=(staf(ordre))';
M=M(:,ordre); % Mmat=Mmat(:,ordre);
ifd=ifd(ordre); ipsm=ipsm(ordre); psm=(psm(ordre))'; 
refhd=refhd(ordre); max_pres=max_pres(ordre);
dparv=(dparv(ordre))';

%Extreme du parametre et niveaux de reference
if strcmp(param,'VGEO'), 
  refhd=min(refhd(1:end-1),refhd(2:end));
  [M2,x2,refhd]=vgeost(pressure(p,mean(lat)),M/100,lat,lg,proj,refhd',max_pres); %%approx!
else M2=M; end;
Mi=min(min(M2(~isnan(M2)))); Ma=max(max(M2(~isnan(M2))));
fprintf('Le parametre choisi varie entre %g et %g .\n',[Mi Ma]);
if strcmp(param,'DYNH') | strcmp(param,'VGEO'),
   fprintf('Les niveaux de ref sont ');
   fprintf('%g\t',refhd); fprintf('\n');
end; %%%

%% Le probleme des bords de la matrice de mesure est traite ci-dessous;
%% des valeurs sont extrapolees pour que les contours se dessinent
%% correctement pres des bords (cad du fond).
%% Prioritairement, l'extrapolation a la station i se fait grace au gradient 
%% horizontal gradh entre la station consecutive la plus profonde (ksta) 
%% et celle 'd'avant' ksta2; ou, si c'est impossible, par repetition du dernier
%% gradient horizontal mesure dparh entre les stations i et ksta.
%% Cette methode peut cependant introduire des abherrations dans le gradient
%% vertical gradvi de la colonne i, soit par comparaison avec le grad vertical
%% gradvk de la colonne ksta (opt 1 et 2), soit par une valeur trop grande
%% (opt 3 et 4). La correction de ces tendances peut se faire chacune de 2
%% facons possibles, ce qui introduit 4 resolutions du probleme (iopt 1,2,3
%% et 4) qui sont toutes les quatre acceptables en principe. Il y aura donc
%% 4 graphes affiches entre lesquelles l'utilisateur pourra choisir.

for ipt=iopt,
 for i=1:ista,
  if ipsm(i)~=lp & ~(strcmp(param,'h')&isinf(ref)),
%%  on trouce ksta, l'indice de la station consecutive la plus profonde
   if i==1, ksta=find(ipsm==max(ipsm(1),ipsm(2)));
            if any(ksta==i), ksta=i; else ksta=2; end;
   elseif i==ista, ksta=find(ipsm==max(ipsm(ista-1),ipsm(ista)));
            if any(ksta==i), ksta=i; else ksta=ista-1; end;
   else     ksta=find(ipsm==max(max(ipsm(i-1),ipsm(i+1)),ipsm(i)));
            if any(ksta==i), ksta=i; 
            else ksta=min(ksta(ksta>=i-1 & ksta<=i+1)); end;
   end;  
   k=ipsm(ksta); kf=ifd(ksta); 
   sig=ksta-i; 		%% sig=-1|0|1 donne le sens de ksta par rapport a i
   ksta2=ksta+sig; 	%% station consecutive a ksta dans le sens de sig
%%%disp(sprintf('%g %g %g %g %g %g %g',[staf(i),sig,i,ksta,ksta2,k,kf]));
%%%disp(sprintf('%g %g %g %g',[ipsm(i),ifd(i),lp,M(ipsm(i),i)]));
   if sig~=0 
    if ksta2>0 & ksta2<=ista,
      k2=ipsm(ksta2);
      rapx=(x(i)-x(ksta))/(x(ksta)-x(ksta2));	%% rapport des distances entre
						%% i, ksta et ksta2
      if k2>=k,
         gradh=(M(ipsm(i)+1:kf,ksta)-M(ipsm(i)+1:kf,ksta2))*rapx;
         M(ipsm(i)+1:kf,i)=M(ipsm(i)+1:kf,ksta)+gradh;
      elseif k2>ipsm(i),
         gradh=(M(ipsm(i)+1:k2,ksta)-M(ipsm(i)+1:k2,ksta2))*rapx;
         M(ipsm(i)+1:k2,i)=M(ipsm(i)+1:k2,ksta)+gradh;
         dparh=M(k2,i)-M(k2,ksta);
         M(k2+1:kf,i)=M(k2+1:kf,ksta)+dparh;
      else
         dparh=M(ipsm(i),i)-M(ipsm(i),ksta);
         M(ipsm(i)+1:kf,i)=M(ipsm(i)+1:kf,ksta)+dparh;
      end;
    else dparh=M(ipsm(i),i)-M(ipsm(i),ksta);
         M(ipsm(i)+1:kf,i)=M(ipsm(i)+1:kf,ksta)+dparh;
    end;
   else if ifd(i)<lp, cs=1; else cs=0; end;
        M(ipsm(i)+1:ifd(i)+cs,i)=...
                  M(ipsm(i),i)...  %*ones(size(M(ipsm(i)+1:ifd(i)+cs,i)));
                  +cumsum(dparv(i)*ones(size(M(ipsm(i)+1:ifd(i)+cs,i))));
%%% disp(sprintf('%g %g %g %g',[ipsm(i),ifd(i)+1,M(ipsm(i),i),M(ifd(i)+1,i)]));
   end;
  elseif strcmp(param,'HDYN')&isinf(ref), 
   M(ipsm(i):lp,i)=zeros(size(M(ipsm(i):lp,i)));
  end;
 end;

 for i=ista:-1:1,
  if ipsm(i)~=lp & ~(strcmp(param,'HDYN')&isinf(ref)),
%%  on trouce ksta, l'indice de la station consecutive la plus profonde
   if i==1, ksta=find(ipsm==max(ipsm(1),ipsm(2)));
            if any(ksta==i), ksta=i; else ksta=2; end;
   elseif i==ista, ksta=find(ipsm==max(ipsm(ista-1),ipsm(ista)));
            if any(ksta==i), ksta=i; else ksta=ista-1; end;
   else     ksta=find(ipsm==max(max(ipsm(i-1),ipsm(i+1)),ipsm(i)));
            if any(ksta==i), ksta=i; 
            else ksta=min(ksta(ksta>=i-1 & ksta<=i+1)); end;
   end;  
   k=ipsm(ksta); kf=ifd(ksta); 
   sig=ksta-i; 		%% sig=-1|0|1 donne le sens de ksta par rapport a i
   ksta2=ksta+sig; 	%% station consecutive a ksta dans le sens de sig
%%%disp(sprintf('%g %g %g %g %g %g %g',[staf(i),sig,i,ksta,ksta2,k,kf]));
%%%disp(sprintf('%g %g %g %g',[ipsm(i),ifd(i),lp,M(ipsm(i),i)]));
   if sig~=0 
    if ksta2>0 & ksta2<=ista,
      k2=ipsm(ksta2);
      rapx=(x(i)-x(ksta))/(x(ksta)-x(ksta2));	%% rapport des distances entre
						%% i, ksta et ksta2
      if k2>=k,
         gradh=(M(ipsm(i)+1:kf,ksta)-M(ipsm(i)+1:kf,ksta2))*rapx;
         M(ipsm(i)+1:kf,i)=M(ipsm(i)+1:kf,ksta)+gradh;
      elseif k2>ipsm(i),
         gradh=(M(ipsm(i)+1:k2,ksta)-M(ipsm(i)+1:k2,ksta2))*rapx;
         M(ipsm(i)+1:k2,i)=M(ipsm(i)+1:k2,ksta)+gradh;
         dparh=M(k2,i)-M(k2,ksta);
         M(k2+1:kf,i)=M(k2+1:kf,ksta)+dparh;
      else
         dparh=M(ipsm(i),i)-M(ipsm(i),ksta);
         M(ipsm(i)+1:kf,i)=M(ipsm(i)+1:kf,ksta)+dparh;
      end;
    else dparh=M(ipsm(i),i)-M(ipsm(i),ksta);
         M(ipsm(i)+1:kf,i)=M(ipsm(i)+1:kf,ksta)+dparh;
    end;
   else if ifd(i)<lp, cs=1; else cs=0; end;
        M(ipsm(i)+1:ifd(i)+cs,i)=...
                  M(ipsm(i),i)*ones(size(M(ipsm(i)+1:ifd(i)+cs,i)));
%                  +cumsum(dparv(i)*ones(size(M(ipsm(i)+1:ifd(i)+cs,i))));
%%% disp(sprintf('%g %g %g %g',[ipsm(i),ifd(i)+1,M(ipsm(i),i),M(ifd(i)+1,i)]));
   end;
   jeg=0;
%%  Comparaison des signes des gradv de la col ksta et de la col i
   for j=ipsm(i)+1:kf,
    gradvk=M(j,ksta)-M(j-1,ksta); 
    gradvi=M(j,i)-M(j-1,i); 
    sig=sign(gradvk);
%% 1-S'ils sont opposes, on force ce gradv au grad de ksta dans col i (opt 1,2)
    if ipt==1 | ipt==2,
     if (sign(gradvk*gradvi)<=0), 
        M(j,i)=M(j-1,i)+gradvk; gradvi=gradvk; jeg=jeg+1; 
     end;
%% 2-On peut aussi forcer ce gradv a 0 (opt 3,4)...
    else
     if (sign(gradvk*gradvi)<=0),
        M(j,i)=M(j-1,i); gradvi=0; jeg=jeg+1; 
     end;
    end;
%% Dans certains cas, il faut empecher le gradv de faire de trop grand sauts ...
    ia=min(j-1-jeg,4);
    gradvia=abs(M(j-jeg-1,i)-M(j-jeg-ia,i))/(ia-1); %moy du gradv au dessus
%% 3-... apres le fond (opt 1,3)
    if ipt==1 | ipt==3,
     if (abs(gradvi)>gradvia) & (abs(gradvi)>abs(gradvk)) & j>ifd(i)...
        & gradvia~=0,
       M(j,i)=M(j-1,i)+sig*gradvia; 
     end;
%% 4-... apres la derniere valeur experimentale (opt 2,4)
    else
     if (abs(gradvi)>gradvia)& (abs(gradvi)>abs(gradvk)) & gradvia~=0,
       M(j,i)=M(j-1,i)+sig*gradvia; 
     end;
    end;
   end;
  end;
 end;

%Calcul si vitesses geostrophiques :
  if strcmp(param,'VGEO'),  %%%
    [M,xpt,refu,refon]=vgeost(pressure(p,mean(lat)),M/100,lat,lg,proj,refhd,max_pres); %%%
    refon=depth(refon,mean(lat));
  else xpt=x;
  end; %%%
 

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mise en forme de la figure
 if ~holdon, figure;clf; set(gcf,'Position',[(20+(ipt-1)*lef) bot larg haut],'COlor','w'); 
 elseif length(iopt)==4, figure(ipt); 
 end;
 if strcmp(ori,'pay'), paysage; else portrait; end;
 pve1=posi(2)/Hfig; dve1=posi(4)/Hfig;
 mha2=pve1-0.04*dve1-0.02; ha2=1.18*dve1+.04;  %+.105;
% mha2=pve1-.2*dve1; ha2=1.5*dve1;  %+.105;
                             
 load mapcolor2;
 ledge=length(plim);
 plim=floor(plim);
 pmin=floor(plim(1)/10)*10; pmax=ceil(plim(ledge)/10)*10;
 dve=dve1*pourc/100;  %hauteur occupee par les futurs axes a1
 pve=fliplr(pve1+cumsum([0 fliplr(dve(2:ledge-1))]));
 if isempty(xlim), xmi=x(1); xma=x(ista);
 else xmi=xlim(1); xma=xlim(2);
 end;
% La figure faisant posi(3) cm en x, 1cm=1/posi(3) unites normalisees
 L=posi(3)/Lfig;               %larg du plot en unite norm.
 if isempty(ech),
      lrg=L;
      ech=posi(3)/(xma-xmi);
      if (proj=='mi' | proj=='km'), ech=ech*100; end;
      fprintf('Echelle : %5.2f cm par unite\n',ech);
 else lrg=ech*(xma-xmi)/Lfig;
      if lrg>L, lrg=L; 
       ech=posi(3)/(xma-xmi);
       if (proj=='mi' | proj=='km'), ech=ech*100; end;
       fprintf('Figure trop grande, echelle ramenee a %5.2f cm par unite\n',ech);
      end;
 end;
 if isfinite(posi(1)), mleft=posi(1)/Lfig; else mleft=(1-lrg)/2; end;

%% Trace (et evaluation des valeurs des contours)
 plim=-plim; p=-p; psm=-psm; lso=-lso;

 for ip=ledge-1:-1:1,
  pmin=plim(ip); pmax=plim(ip+1);
%  M2=M(p>=pmin & p<=pmax+dp,:); x2=xpt; p2=p(p>=pmin & p<=pmax+dp);
  M2=M(p<=pmin & p>=pmax-dp,:); x2=xpt; p2=p(p<=pmin & p>=pmax-dp);
  M2(abs(M2)<1e-4)=zeros(size(M2(abs(M2)<1e-4)));
  lx=length(x2);
  M2i=min(min(M2(~isnan(M2)))); M2a=max(max(M2(~isnan(M2))));

%% Mise en forme des axes
  a1(ip)=axes('Position',[mleft pve(ip) lrg dve(ip)]);
%  set(a1(ip),'YDir','reverse','YLim',[pmin pmax]);
  set(a1(ip),'YLim',[pmax pmin]);
  set(a1(ip),'FontSize',10,'XLim',[xmi xma],'Box','on'); grid on;
  li=line([xmi;xmi;xma;xma;xmi],[pmin;pmax;pmax;pmin;pmin]);
  set(li,'Color','k');
  hold on
%% Mise en couleur
  if col==1,
    %x1=divrond(x2,100); 
    %p1=divrond(p2,100); 
    %M1=griddata(x2',p2,M2,x1',p1);
    x1=x2;
    p1=p2;
    M1=M2;
    su=pcolor(x1',p1,M1); set(su,'ZData',-1+0.*M1);
    shading flat;
    caxis([Mi Ma]); colormap(jetpers);
  end;
  %x1=divrond(x2,1000); 
  %M2=griddata(x2',p2,M2,x1',p2,'cubic');
  %x2=x1;
  if ~isempty(fil)&length(fil(:))>2, continterfill(x2,p2,M2,fil(:,1:2),fil(:,3:5)); end;
  if strcmp(cont,'on'),
%  if ~isempty(dpa) | ~isempty(vpa),
    if isempty(dpa) & isempty(vpa), vpa=divrond([Mi Ma],10); disp(vpa'); end;
    if ~isempty(dpa),
     odg=10.^(round(log10(dpa(ip))));
     pami=floor(M2i./odg).*odg;
     va=(pami:dpa(ip):M2a)';
    else va=vpa(~isnan(vpa(:,ip)),ip);
    end;
%%  Trace des contours gras
    if ~isempty(vpb),
     vb=vpb(~isnan(vpb(:,ip)),ip);
     if length(vb)==1, vb=[vb vb]; end;
     if any(vb>=M2i & vb<=M2a), 
      [Ccb,hcb]=contour(x2,p2,M2,vb','k-'); set(hcb,'LineWidth',2);
      if strcmp(contlab,'on'),
        labcb=clabel(Ccb,hcb); set(labcb,'FontSize',7,'Color','k','FontWeight','bold');
      end
     else fprintf('Aucun contour gras dans la %geme partie\n',ip);
     end;
%%   On enleve les valeurs = a vpb dans vpa
     va2=sort([va;vb]);
     ivpa=find(diff(va2)==0);ivpa=sort([ivpa;ivpa+1]);
     va2(ivpa)=[];
     va=sort([va;va2]); ivpa=find(diff(va)==0);
     va=va(ivpa);
     if length(va)==1, va=[va va]; end;
    end;
%%  Trace des contours fins
    if any(va>=M2i & va<=M2a),
     if ~isempty(fil)&length(fil(:))==2,
      [Ccaf,hcaf]=contourf(x2,p2,M2,va'); shading flat;
      [Cca,hca]=contour(x2,p2,M2,va','k-');
%      labca=clabel(Cca,hca); set(labca,'FontSize',7,'Color','k','FontWeight','bold');
      if strcmp(contlab,'on'),
        labca=clabel(Cca,hca,'FontSize',7,'Color','k','FontWeight','bold','LabelSpacing',300);
      end
      if isnan(fil(1)), fil(1)=Mi; end; if isnan(fil(2)), fil(2)=Ma; end;
      colormap(clr1); caxis(fil);
     else
      [Cca,hca]=contour(x2,p2,M2,va','k-');
      Ccaf=[];hcaf=[];
      if strcmp(contlab,'on'),
        labca=clabel(Cca,hca); set(labca,'FontSize',7,'Color','k','FontWeight','bold');
      end
     end;
    else fprintf('Aucun contour fin dans la %geme partie\n',ip);
    end;
%  else
%   [Cc,hc]=contour(x2,p2,M2);
%   labc=clabel(Cc,hc); set(labc,'Color','k','FontSize',7,'FontWeight','bold');
%  end;
  end;
%  if ~isempty(fil)&length(fil(:))>2, continterfill(x2,p2,M2,fil(:,1:2),fil(:,3:5)); end;
%% trace du cache blanc pour les vitesses geostrophiques
  if strcmp(param,'VGEO'),
   refon=-refon; 
   [icach,ycach]=stairs([refon(:);refon(lx)]);
%   patch([xmi;x(icach+1);xma],[pmax;ycach;pmax],[1 1 1],'LineWidth',.5);
   patch([xmi;x(icach);xma],[pmax;ycach;pmax],[1 1 1],'LineWidth',.5);
  end;
%%efface sous les dernieres mesures et trace la ligne des dernieres mesures:
  pmes=patch([xmi;xmi;x;xma;xma],[pmax;max(psm);psm;max(psm);pmax],[1 1 1],'EdgeColor','none'); 
  set(pmes,'ZData',get(pmes,'XData')*0+9);
  lmes=line([xmi;x;xma],[max(psm);psm;max(psm)],[xmi;x;xma]*0+9,'LineStyle','--'); 
%% trace du fond en gris
  pfond=patch([xmi;x;xma],[pmax;lso;pmax],[0.8 0.8 0.8],'LineWidth',.5);
  set(pfond,'ZData',get(pfond,'XData')*0+10);
%  line([xmi;x;xma],[pmax;lso;pmax]);

%% Ecriture des valeurs de la matrice sur la coupe
  if strcmp(mat,'on'),
   p3=(p(min(find(p>=plim(ip)))):dmat:p(max(find(p<=plim(ip+1)))))';
   for i=1:lx,    
    if dmat==dp, M3=M2(:,i); p3=p2; else M3=interp1(p2,M2(:,i),p3); end;
    lM3=length(M3); M3lab=[];
    for j=1:lM3, M3lab=str2mat(M3lab,sprintf('%8.4f',M3(j))); end;
    M3lab=M3lab(2:lM3+1,:);
    tM3=text(x(i)*ones(lM3,1),p3,M3lab,'FontSize',5);
    set(tM3,'Horizo','center','vertical','middle');
    line([xmi;x;xma],[max(psm);psm;max(psm)]);
   end;
  end;

%% Mise en forme du dessin
  if ip~=ledge-1, set(a1(ip),'XTickLabel',[]); end;
%% Ticks correspondant aux stations
  tista=line([x';x'],[pmin;pmin+(pmax-pmin)/100]*ones(size(x')),'Color','k');

%% reecriture des labels de pression en positif
  yti=divrond([pmax pmin],8);
  ylab=[]; for ila=1:length(yti), ylab=strvcat(ylab,num2str(-yti(ila))); end;
  set(gca,'YTick',yti,'YTickLabel',ylab);
 end;

%% Ecriture des stations au dessus du graphe
 xlab=[];
 for i=1:ista, xlab=str2mat(xlab,sprintf('%g',staf(i))); end;
 blk=xlab(1,:);
 xlab=xlab(2:ista+1,:);
 if     ista>40, for i=2:2:ista, xlab(i,:)=blk; end;
 elseif ista>80, for i=3:4:ista, xlab(i,:)=blk; end;
 end;
 testa=text(x,pmin*ones(size(x)),xlab,'Horizo','center','Color','k');
 set(testa,'Vertical','bottom','FontSize',8,'FontAngle','italic');

 axes(a1(ledge-1));
 xlabel(xlabe);
 ylabel('Depth (m)');
%% Textes
% a2=axes('Position',[mleft 0.1 lrg 0.8],'Visible','off');
 a2=axes('Position',[mleft mha2 lrg ha2],'Visible','off');
 if ~isempty(titre) | ~isempty(com),
   if ~isempty(titre),
    ttitre=text(0.5,.97,titre,'Horizo','center','FontWeight','bold','FontSize',12);
   end
   if ~isempty(com),
    tcom=text(0.5,0.91,com,'Horizo','center','FontSize',10);
   end
 end;

 if larg==250,
  topt=text(0.5,0,['option ',int2str(ipt)],'Horizo','center','FontSize',8);
 end;

 tdate=text(1,0.89,dattext,'Horizo','right','FontSize',8);

 if strcmp(param,'h') | strcmp(param,'u'),
    if isinf(ref), reftext='reference = plus grande pression atteinte';
    else reftext=['Reference = ',num2str(ref),' dbar'];
    end;
    tref=text(0,0.88,reftext,'Horizo','left','FontSize',8);
 end;
 
 if islg==1,     tNSEW=text([0;0.98],[0.02;0.02],['W';'E'],'FontWeight','bold');
 elseif islg==0, tNSEW=text([0;0.98],[0.02;0.02],['S';'N'],'FontWeight','bold');
 end;
 axes(a1(ledge-1));
end;
%set(gcf,'PaperPosition',[1 1 27 19]);

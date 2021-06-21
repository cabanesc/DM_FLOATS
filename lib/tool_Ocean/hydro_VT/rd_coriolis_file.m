
function [psal,temp,tpot,sig0,vorp,pres,vorpi,juld,platform_number,cycle_number,latitude,longitude]=rd_coriolis_file(inpath,fname,std_AREA,typ_dta)

% type_dta='rt' ou 'dm'

ipsal = 0;ipres =0;ideph=0;itemp=0;

%ncstartup
f_cdf = netcdf([inpath fname],'read');

% Reads general informations

platform_number  = f_cdf{'PLATFORM_NUMBER'}(:);
cycnum=f_cdf{'CYCLE_NUMBER'}(:);

latitude= f_cdf{'LATITUDE'}(:);
longitude= f_cdf{'LONGITUDE'}(:);
position_qc  = f_cdf{'POSITION_QC'}(:);
juld=f_cdf{'JULD'}(:);
juld_qc      = f_cdf{'JULD_QC'}(:);
n_prof       = length(latitude);


% Selection des donnees sur critere de zone geographique et de QC
iok_geog = find(latitude  > std_AREA(1) & latitude  < std_AREA(2) & ...
                longitude > std_AREA(3) & longitude < std_AREA(4) );
%fprintf(std_mess,[' Nombre de profils lus : %d\n',...
%       ' Nombre de profils dans la zone fixee (iok_geog): %d\n'],...
%       n_prof, length(iok_geog));

iok_qc_pos = find(double(position_qc)==double('0') | ...
                  double(position_qc)==double('1') | ...
                  double(position_qc)==double('2') | ...
                  double(position_qc)==double('5'));
iok_juld_qc = find(double(juld_qc)==double('0') | ...
                   double(juld_qc)==double('1') | ...
                   double(juld_qc)==double('2') | ...
                   double(juld_qc)==double('5'));
iok_qc_pos_date = intersect(iok_qc_pos,iok_juld_qc);
iok10 = intersect(iok_geog,iok_qc_pos_date);

% Boucle 10 
if  ~isempty(iok10)

  latitude      = latitude(iok10);
  longitude     = longitude(iok10);
  data_type     = f_cdf{'DATA_TYPE'}(:);
  format_version= f_cdf{'FORMAT_VERSION'}(:);
  reference_date_time = f_cdf{'REFERENCE_DATE_TIME'}(:);
  bid = f_cdf{'PROJECT_NAME'}(:);project_name = bid(iok10,:);
  bid = f_cdf{'PI_NAME'}(:);pi_name = bid(iok10,:);
  bid = f_cdf{'PLATFORM_NUMBER'}(:);platform_number = bid(iok10,:);
  bid = f_cdf{'CYCLE_NUMBER'}(:);cycle_number= bid(iok10);
  bid = f_cdf{'DIRECTION'}(:);direction  =bid(iok10);
  bid = f_cdf{'DATA_CENTRE'}(:);data_centre = bid(iok10,:);
  bid = f_cdf{'DATA_STATE_INDICATOR'}(:);data_state_indicator =bid(iok10,:);
  bid = f_cdf{'DATA_MODE'}(:);data_mode = bid(iok10);
  bid = f_cdf{'DC_REFERENCE'}(:);dc_reference = bid(iok10,:);
  bid = f_cdf{'INST_REFERENCE'}(:);inst_reference = bid(iok10,:);
  bid = f_cdf{'WMO_INST_TYPE'}(:);wmo_inst_type = bid(iok10,:);
  bid = f_cdf{'JULD'}(:);juld = bid(iok10);
  
  clear bid
  
  % Lecture des mesures, selection en fonction des QCs
  
  tabparam=['PRES';'DEPH';'TEMP';'PSAL'];
  if strcmp(typ_dta,'rt') == 1
    tabparam2=tabparam;
  elseif strcmp(typ_dta,'dm') == 1
    tabparam2=['PRES_ADJUSTED';'DEPH_ADJUSTED';'TEMP_ADJUSTED';'PSAL_ADJUSTED'];
  end


  for iparam =1:4
    vparam=tabparam(iparam,:);
    vparam2=tabparam2(iparam,:);
    bid = f_cdf{[vparam]}(:);
    if ~isempty(bid)
      disp(['traitement de ' vparam2])
      param=bid(iok10,:);
      bid = f_cdf{[vparam2 '_QC']}(:);param_ini_qc=bid(iok10,:);
      
      if strcmp(vparam2,'DEPH') == 1  % pas de FillValue dans les fichiers !
        inok = find(param<0 | param>9900);
      else
        param_fillval = f_cdf{[vparam2]}.FillValue_(:);
        inok = find(param==param_fillval);
      end
      if ~isempty(inok)
        param(inok)=NaN;
      end
      if strcmp(vparam2,'DEPH') == 1 | strcmp(vparam2,'PRES') == 1 | strcmp(vparam2,'PRES_ADJUSTED') == 1
        inok =find(double(param_ini_qc)~= double('0') &...
                   double(param_ini_qc)~= double('1') &...
                   double(param_ini_qc)~= double('2') );
      else
             
             
        inok =find(double(param_ini_qc)~= double('1') &...
                   double(param_ini_qc)~= double('2'));
      end     
      if ~isempty(inok)
        param(inok)=NaN;
      end
      iok=find(isfinite(param));
      eval([lower(vparam) '=param;']);
      [nx,ny]=size(param_ini_qc);
      for ix=1:nx
        for iy=1:ny
          if strcmp(param_ini_qc(ix,iy),' ') == 1
            qcnum(ix,iy)=NaN;
          else
            qcnum(ix,iy)=sscanf(param_ini_qc(ix,iy),'%f');
          end
        end
      end
      eval([lower(vparam) '_ini_qc=qcnum;']);
      if ~isempty(iok)
        eval(['i' lower(vparam) '=1'])
      else
        %fprintf(std_mess,'%s\n',...
        %  ['  Aucune temperature valide dans ce fichier']);
      end
    end
    clear param_ini_qc qcnum param vparam
  end
  

  % On ne continue qu'a condition d'avoir TEMP et PSAL et PRES ou DEPH
  % (pour calculer sig0 apres)  
  % Boucle 15
  if (itemp==1 & ipsal==1) & (ipres == 1 | ideph == 1)
    
    n_prof=length(iok10);
    
    % si une mesure de temp, psal, pres ou deph =NaN:
    % les autres parametres correspondants sont mis a NaN
     
    inok1=find(~isfinite(temp) | ~isfinite(psal));
    if ipres==1 && ideph==1
      inok2=find(~isfinite(pres) & ~isfinite(deph));
    elseif ipres==1 && ideph==0
      inok2=find(~isfinite(pres));
    elseif ipres==0 && ideph==1
      inok2=find(~isfinite(deph));
    end
    inok=union(inok1,inok2);
    if ~isempty(inok)
      if(ipres==1)
        pres(inok)=NaN;
      end
      if(ideph==1)
        deph(inok)=NaN;
      end
      psal(inok)=NaN;
      temp(inok)=NaN;
    end

    %nvalok_temp: nb de mes. valides pour chaque profil en temp
    nvalok_temp = NaN*ones(n_prof,1);
    for i=1:n_prof % 16
      t=squeeze(temp(i,:));
      iok=find(isfinite(t));
      if ~isempty(iok)
        nvalok_temp(i) = length(iok);
      else
        nvalok_temp(i) = 0;
      end
    end %16
    iok=[];
  
    % On ne garde pas les profils qui ont aucune mes. valide
    % en temp
    iok11 = find(nvalok_temp>0);
    if ((~isempty(iok11)) & (length(iok11)~=n_prof)) % 18
      latitude=latitude(iok11);
      longitude=longitude(iok11);
      project_name=project_name(iok11,:);
      pi_name=pi_name(iok11,:);
      platform_number=platform_number(iok11,:);
      cycle_number=cycle_number(iok11);
      direction=direction(iok11);
      data_centre=data_centre(iok11,:);
      data_state_indicator=data_state_indicator(iok11,:);
      data_mode=data_mode(iok11);
      dc_reference=dc_reference(iok11,:);
      inst_reference=inst_reference(iok11,:);
      wmo_inst_type=wmo_inst_type(iok11,:);
      juld=juld(iok11);
      temp=temp(iok11,:);
      psal=psal(iok11,:);
      temp_ini_qc=temp_ini_qc(iok11,:);
      psal_ini_qc=psal_ini_qc(iok11,:);
      nvalok_temp=nvalok_temp(iok11);
      if ipres==1
        pres=pres(iok11,:);
        pres_ini_qc=pres_ini_qc(iok11,:);
      end
      if ideph==1
        deph=deph(iok11,:);
        deph_ini_qc=deph_ini_qc(iok11,:);
      end
      n_prof=length(nvalok_temp);
    end % fin 18
  else
    n_prof=0; 
  end 
  % Fin boucle 15, test sur itemp et ipsal
     % s'il y a des profils, alors on continue
  % Boucle 19
  if n_prof>0
    
    pres_new=NaN*ones(size(temp));
    
    % On boucle sur les profils
    % Boucle 20    
    for i_prof=1:n_prof
      iok=[];
      ideph_prof=0;ipres_prof=0;
      ipsal_prof=0;itemp_prof=0;
      
      % On verifie qu'il y a des donnees de pression ou d'immersion
      % on considere que c'est ok
     
      t=squeeze(temp(i_prof,:));
      s=squeeze(psal(i_prof,:));
      dorp=squeeze(pres(i_prof,:));
      tqc=squeeze(temp_ini_qc(i_prof,:));
      sqc=squeeze(psal_ini_qc(i_prof,:));
      pqc=squeeze(pres_ini_qc(i_prof,:));
      lat=latitude(i_prof);
      
      
      % On enleve les NaN sur t et dorp
      iokt=find(isfinite(t)&isfinite(dorp));
      dorp=dorp(iokt);
      t=t(iokt);
      s=s(iokt);
      tqc=tqc(iokt);
      sqc=sqc(iokt);
      pqc=pqc(iokt);
      ioks=find(~isfinite(s));
      if isempty(ioks) == 0
       disp(['PSAL = NaN ']);
       break;
      end  
      iokd=find(~isfinite(dorp));
      if isempty(iokd) == 0
        disp(['PRES = NaN ']);
        break;
      end
      
      nvalok_temp(i_prof)=length(t);
      
      %110 ++++++++++++++++++++++++++++++++++++++++++
      pres_new(i_prof,1:length(dorp))=dorp;
      temp(i_prof,:)=NaN;
      temp(i_prof,1:length(t))=t;
      psal(i_prof,:)=NaN;
      psal(i_prof,1:length(t))=s;
      temp_ini_qc(i_prof,:)=NaN;
      temp_ini_qc(i_prof,1:length(t))=tqc;
      psal_ini_qc(i_prof,:)=NaN;
      psal_ini_qc(i_prof,1:length(t))=sqc;
      pres_ini_qc(i_prof,:)=NaN;
      pres_ini_qc(i_prof,1:length(t))=pqc;
    end
    % Fin boucle 20 for i_prof=1:n_prof
 
    %suppression des profils entierement a NaN
    %suppression des profils entierement a NaN
    % par ex ceux pour lesquels il n y avait
    % aucun sig0.
    iok12=find(nvalok_temp>0);
    if ((~isempty(iok12)) & (length(iok12)~=n_prof))
      latitude=latitude(iok12);
      longitude=longitude(iok12);
      project_name=project_name(iok12,:);
      pi_name=pi_name(iok12,:);
      platform_number=platform_number(iok12,:);
      cycle_number=cycle_number(iok12);
      direction=direction(iok12);
      data_centre=data_centre(iok12,:);
      data_state_indicator=data_state_indicator(iok12,:);
      data_mode=data_mode(iok12);
      dc_reference=dc_reference(iok12,:);
      inst_reference=inst_reference(iok12,:);
      wmo_inst_type=wmo_inst_type(iok12,:);
      juld=juld(iok12);
      temp=temp(iok12,:);
      nvalok_temp=nvalok_temp(iok12);
      pres_new=pres_new(iok12,:);
      psal=psal(iok12,:);
      pres_ini_qc=pres_ini_qc(iok12,:);
      temp_ini_qc=temp_ini_qc(iok12,:);
      psal_ini_qc=psal_ini_qc(iok12,:);
      n_prof=length(nvalok_temp);
    elseif isempty(iok12)
      n_prof=0;
    end
  end %if n_prof>0 5 boucle 19

pres=pres_new;
tpot=tetai(pres,temp,psal,0);
[svan,sig0]=swstat90(psal,tpot,0);
size(psal)

ns=size(sig0);
vorpi=NaN*ones(ns);
for i=1:ns(1)
  frequence_B_V=fbruva(psal(i,:)',temp(i,:)',pres(i,:)',latitude(i));
  vrp=vorpot(frequence_B_V,latitude(i));
  vorp(i,:)=vrp';

  presi=[0:10:2000];
  ifin=find(isfinite(pres(i,:)));
  if length(ifin)>=2  
    tempi=interp1(pres(i,ifin)',temp(i,ifin)',presi');
    psali=interp1(pres(i,ifin)',psal(i,ifin)',presi');
    templ=lanczos(tempi',1/20,11);
    psall=lanczos(psali',1/20,11);
    fbvi=fbruva(psall,templ,presi',latitude(i));
    vrpi=vorpot(fbvi,latitude(i));
    %vorpll=lanczos(vrpi',1/20,11);
    %vrpi=vorpll;
    ifin2=find(isfinite(vrpi));
    if length(ifin2)>=2   
      vrpf=interp1(presi(ifin2)',vrpi(ifin2),pres(i,ifin)');
      vorpi(i,ifin)=vrpf';
    end     
  end
end

else
pres=[];
tpot=[];
vorp=[];
temp=[];
psal=[];
sig0=[];
end


close(f_cdf)



function [cnews] = routine_profil1_cpcor(fw,nocycl,floatnum,po_system_configuration,palier,pas_surface,pas_fond)

%==========================================================================
% OBJECTIF: comparer un profil Argo a un profil CTD de reference (shipboard
% reference CTD), et calculer le coefficient cpcor optimun et
% éventuellement un factueur multiplicatif pour la conductivité
% (conservative conductivity)
%  programme basé sur le programme Co_press_corr.m de Greg Johnson.
%==========================================================================
global CPcor
global CTcor
global temp_argo pres_argo psal_argo lon_argo lat_argo psal_campagne ct_campagne pres_campagne
global minDEPTH
global press_offset
global CPcor_fixed
global M0




% nominal CT values from SBE.
%--------------------------------------------------------------------------
CPcor = -9.5700E-8;
CTcor =  3.2500E-6;

% Minimum depth considered for the optimum fit
minDEPTH = po_system_configuration.ZOOM

% starting values for optimisation function
inp=[-13.0e-8  1 ];
%inp=[-0.1715  0];
options = optimset('TolFun',1e-9);

%--------------------------------------------------------------------------
% lecture du profil Argo spécifié
%--------------------------------------------------------------------------

repnc=[po_system_configuration.FLOAT_SOURCE_NETCDF po_system_configuration.DACNAME '/' floatnum '/profiles/'];

% lit les 10 premiers cycles
%---------------------------
if isempty(po_system_configuration.DACNAME)==0
    repnc = [po_system_configuration.FLOAT_SOURCE_NETCDF po_system_configuration.DACNAME '/' floatnum '/profiles/'];
else
    repnc = [po_system_configuration.FLOAT_SOURCE_NETCDF po_system_configuration.DACNAME '/' floatnum '/'];
end
if isempty(po_system_configuration.DACNAME)&&isempty(floatnum)
    repnc = [po_system_configuration.FLOAT_SOURCE_NETCDF];calseries
end
IncludeDescProf=1;
[file_list]=select_float_files_on_ftp(floatnum,po_system_configuration.DACNAME,po_system_configuration.FLOAT_SOURCE_NETCDF,'C',IncludeDescProf);
if length(file_list) >= str2double(po_system_configuration.MAX_PROF)
    file_list=file_list(1:str2double(po_system_configuration.MAX_PROF));
end

[F,Dim,file_list] = create_multi_from_filelist(floatnum,po_system_configuration.DACNAME,po_system_configuration.FLOAT_SOURCE_NETCDF,file_list,'Primary',[]);

nocycl=str2num(nocycl);

%--------------------------------------------------------------------------
% lecture du fichier technique (récupère la presssion de surface du cycle1
%--------------------------------------------------------------------------

FILENAME_TECH = [po_system_configuration.FLOAT_SOURCE_NETCDF po_system_configuration.DACNAME '/' floatnum '/' floatnum '_tech.nc'];
T = read_netcdf_allthefile(FILENAME_TECH);
% trouve les pressions de surface
isurf=find(findstr_tab(cellstr(T.technical_parameter_name.data),'PRES_SurfaceOffsetCorrectedNotResetNegative_1cBarResolution_dbar')==1);
cycl=T.cycle_number.data(isurf);
cyc1=cycl==1;
if nocycl==1|nocycl==2
press_offset= str2num(T.technical_parameter_value.data(isurf(cyc1),:));
else
press_offset=0;
end
press_offset
if nocycl<=length(F.cycle_number.data)
    %close all
    F=replace_fill_bynan(F);
    F=format_flags_char2num(F);
    
    nocycl_str=[num2str(F.cycle_number.data(nocycl)) strtrim(F.direction.data(nocycl))];
    lat_argo = double(F.latitude.data(nocycl,:));
    lon_argo = double(F.longitude.data(nocycl,:));
    juld_argo = double(F.juld.data(nocycl,:));
    
    psal_argo = double(F.psal.data(nocycl,:));
    temp_argo = double(F.temp.data(nocycl,:));
    pres_argo = double(F.pres.data(nocycl,:));
    
    psalqc_argo = F.psal_qc.data(nocycl,:)';
    tempqc_argo = F.temp_qc.data(nocycl,:)';
    presqc_argo = F.pres_qc.data(nocycl,:)';
    
    psal_argo(psalqc_argo>3)=NaN;
    
    % compute derived quantities including absolute salinity, sa, conservative
    % temperature, ct, and conservative conductivity, cco.
    %--------------------------------------------------------------------------
    sa_argo = gsw_SA_from_SP(psal_argo,pres_argo,lon_argo,lat_argo);
    ct_argo = gsw_CT_from_t(sa_argo,temp_argo,pres_argo);
    cco_argo = gsw_C_from_SP(psal_argo,ct_argo,0*pres_argo);
    
    % F.ptmp.data = sw_ptmp(F.psal.data,F.temp.data,F.pres.data,0)
    
    nocyclstr = [num2str(F.cycle_number.data(nocycl)) strtrim(F.direction.data(nocycl))];
    
    
    %--------------------------------------------------------------------------
    % Recherche des donnees du profil Campagne le plus proche (geographiquement)
    % du profil
    %--------------------------------------------------------------------------
    
    
    filename = fullfile( po_system_configuration.DATA_DIRECTORY, strcat(po_system_configuration.CAMPAGNE_MAT, po_system_configuration.POSTFIX ))
    
    if strcmp(po_system_configuration.CAMPAGNE_MAT,'PedroVelez')
        filename = fullfile( po_system_configuration.DATA_DIRECTORY, [floatnum ,'.mat']);
    end
    
    lo_float_source_data = load( filename );

    if strcmp(po_system_configuration.CAMPAGNE_MAT,'PedroVelez')
        iok=1;
        profil_ok = 1;
        lat_campagne = double(lo_float_source_data.LATI(profil_ok));
        lon_campagne = double(lo_float_source_data.LONG(profil_ok));
        %juld_campagne = lo_float_source_data.gtime(lo_float_source_data.gtime);
        psal_campagne = lo_float_source_data.salt';
        temp_campagne = lo_float_source_data.temp';
        pres_campagne = lo_float_source_data.pres';
    else
        lat_campagne = lo_float_source_data.LAT;
        lon_campagne = lo_float_source_data.LONG;
        
        dista1 = andoyer(lon_campagne,lat_campagne, lon_argo,lat_argo);
        [distascend,isort] = sort(dista1);
        iok=1;
        profil_ok = isort(iok);
        lat_campagne = double(lo_float_source_data.LAT(profil_ok));
        lon_campagne = double(lo_float_source_data.LONG(profil_ok));
        %juld_campagne = lo_float_source_data.DATES(profil_ok);
        psal_campagne = lo_float_source_data.SAL(profil_ok,:);
        temp_campagne = lo_float_source_data.TEMP(profil_ok,:);
        pres_campagne = lo_float_source_data.PRES(profil_ok,:);
    end
    % compute derived quantities for the reference station including
    % conservative temperature and conservative conductivy (cco1), which is
    % conductivity evaluated using the conservative temperature and assuming
    % zero pressure.
    
    sa_campagne = gsw_SA_from_SP(psal_campagne,pres_campagne,lon_campagne,lat_campagne);
    ct_campagne = gsw_CT_from_t(sa_campagne,temp_campagne,pres_campagne);
    cco_campagne = gsw_C_from_SP(psal_campagne,ct_campagne,0*pres_campagne);
    psal_campagne_r = gsw_SP_from_C(cco_campagne,ct_campagne,0*pres_campagne); %test
    
    % interpolation on float theta's levels
    [psal_i_campagne,pres_i_campagne] = interp_climatology(psal_campagne',ct_campagne',pres_campagne',psal_argo,ct_argo,pres_argo); % routine OW
    cco_i_campagne = gsw_C_from_SP(psal_i_campagne,ct_argo',0*pres_i_campagne);
    difference_pres_theta = [pres_argo-pres_i_campagne'];
    ip = find(difference_pres_theta<150);
    difference_psal_theta = psal_argo-psal_i_campagne';
    
    % on cherche une valeur optimale pour le cpcor et l'offset => resultat
    % dans cnew.
	
	 % on cherche une valeur optimale pour l'offset en utilisant la valeur
    % cpcor nominale
    CPcor_fixed = CPcor;
    inp2=[1];M0=1;
    cnew2 = fminsearch(@myfun,inp2,options)
	%psal_argo_nominal = change_cpcor([CPcor,cnew2(1)]);
	
	% on utilise la valeur de cet offset pour faire une première correction 
    M0=cnew2(1)
    cnew = fminsearch(@myfun,inp,options)
	keyboard
    
    % on cherche une valeur optimale pour l'offset en utilisant la valeur
    % cpcor nominale
    CPcor_fixed = CPcor;
    inp2=[1];
    cnew2 = fminsearch(@myfun,inp2,options)
	
    % on cherche une valeur optimale pour l'offset en utilisant une valeur
    % mediane pour CPCOR=-1.294e-7
    CPcor_fixed = -1.41e-7;
    inp3=[1];
    cnew3 = fminsearch(@myfun,inp3,options)
    
    
    pres_argo_corr=pres_argo + press_offset;
    
   difference_psal_theta_raw = [psal_argo-psal_i_campagne'];
   
   cond_argo = gsw_C_from_SP(psal_argo,temp_argo,pres_argo);
   psal_argo_corr = gsw_SP_from_C(cond_argo,temp_argo,pres_argo_corr);
    
   difference_psal_theta_corr = [psal_argo_corr-psal_i_campagne'];
    

    psal_argo_nominal = change_cpcor([CPcor,cnew2(1)]);
    sa_argo = gsw_SA_from_SP(psal_argo_nominal,pres_argo_corr,lon_argo,lat_argo);
    ct_argo2 = gsw_CT_from_t(sa_argo,temp_argo,pres_argo_corr);
    [psal_i_campagne,pres_i_campagne]=interp_climatology(psal_campagne',ct_campagne',pres_campagne',psal_argo_nominal,ct_argo2,pres_argo_corr); % routine OW
    difference_psal_theta_nominal=[psal_argo_nominal-psal_i_campagne'];
    
    psal_argo_optim = change_cpcor([cnew(1),cnew(2)]);
    sa_argo = gsw_SA_from_SP(psal_argo_optim,pres_argo_corr,lon_argo,lat_argo);
    ct_argo2 = gsw_CT_from_t(sa_argo,temp_argo,pres_argo_corr);
    [psal_i_campagne_o,pres_i_campagne_o]=interp_climatology(psal_campagne',ct_campagne',pres_campagne',psal_argo_optim,ct_argo2,pres_argo_corr); % routine OW
    difference_psal_theta_optim=[psal_argo_optim-psal_i_campagne_o'];
    
    psal_argo_median = change_cpcor([CPcor_fixed,cnew3(1)]);
    sa_argo = gsw_SA_from_SP(psal_argo_median,pres_argo_corr,lon_argo,lat_argo);
    ct_argo2 = gsw_CT_from_t(sa_argo,temp_argo,pres_argo_corr);
    [psal_i_campagne_m,pres_i_campagne_m] = interp_climatology(psal_campagne',ct_campagne',pres_campagne',psal_argo_median,ct_argo2,pres_argo_corr); % routine OW
    difference_psal_theta_median = [psal_argo_median-psal_i_campagne_m'];
    
    
    
    theoffsset_pres = meanoutnan(psal_argo-psal_argo_corr);

    theoffsset_nom = meanoutnan(psal_argo_corr-psal_argo_nominal);
    theoffsset_op = meanoutnan(psal_argo_corr-psal_argo_optim);
    theoffsset_med = meanoutnan(psal_argo_corr-psal_argo_median);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dir_enregistre=[po_system_configuration.PROF1_DIRECTORY '/' floatnum ];
    if exist(dir_enregistre)==0
        mkdir(dir_enregistre)
    end
    
	PRINT=0;
	if PRINT==1
    figure(1)
    hold on
    grid on
    plot(difference_psal_theta_corr,pres_argo_corr,'ok')
    plot(difference_psal_theta_nominal,pres_argo_corr,'+b')
    plot(difference_psal_theta_optim,pres_argo_corr,'+r')
    
    maxval=max(abs(difference_psal_theta_raw(pres_argo>str2num(minDEPTH))));
    if isempty(maxval)||isnan(maxval)
         maxval=max(abs(difference_psal_theta_median(pres_argo>str2num('1000'))));
     end
	if isempty(maxval)||isnan(maxval)
         maxval=max(abs(difference_psal_theta_median(pres_argo>str2num('500'))));
     end
    set(gca,'XLim',[-maxval-0.001,maxval+0.001])
    set(gca,'Ydir','reverse')
    title({['Salinity deviation from the shipboard reference. Float ' floatnum ', Cycle ' nocyclstr ]},'Fontsize',14)
    xlabel('Salinity deviation on float theta levels','Fontsize',12)
    legend({'Original  profile (press corrected): nominal CPCOR  (-9.57e-8) , no offset corrected',['Modified profile: nominal CPCOR  (-9.57e-8) , optimal offset  (' num2str(theoffsset_nom) ')' ],['Modified profile: optimal CPCOR (', num2str(cnew(1)), ') , optimal offset (', num2str(theoffsset_op), ')']},'location','SouthOutside','Fontsize',14,'FontWeight','bold')
    ylabel('pressure','Fontsize',12)
    set(gcf,'papertype','A4','paperunits','centimeters','paperorientation','portrait','paperposition',[.25,.75,20,25]);
    eval(['print -dpdf ' dir_enregistre '/CPCOR_analysis_' floatnum '_' nocyclstr   '_new.pdf']);
    
    figure(2)
    hold on
    grid on
    plot(difference_psal_theta_raw,pres_argo,'.k')
    plot(difference_psal_theta_corr,pres_argo_corr,'ok')
    set(gca,'Ydir','reverse')
    
    maxval=max(abs(difference_psal_theta_raw(pres_argo>str2num(minDEPTH))));
    if isempty(maxval)||isnan(maxval)
         maxval=max(abs(difference_psal_theta_median(pres_argo>str2num('1000'))));
     end
	if isempty(maxval)||isnan(maxval)
         maxval=max(abs(difference_psal_theta_median(pres_argo>str2num('500'))));
     end
    set(gca,'XLim',[-maxval-0.001,maxval+0.001])
    
    title({['Salinity deviation from the shipboard reference. Float ' floatnum ', Cycle ' nocyclstr ]},'Fontsize',14)
    xlabel('Salinity deviation on float theta levels','Fontsize',12)
    legend({'Original  profile: nominal CPCOR  (-9.57e-8) , no offset corrected',['Original  profile with pressure offset corrected (' num2str(press_offset) 'db), salinity is re-calculated']},'location','SouthOutside','Fontsize',14,'FontWeight','bold')
    ylabel('pressure','Fontsize',12)
    set(gcf,'papertype','A4','paperunits','centimeters','paperorientation','portrait','paperposition',[.25,.75,20,25]);
    eval(['print -dpdf ' dir_enregistre '/press_correction_cycle1_' floatnum '_' nocyclstr   '_new.pdf']);
    
    figure(3)
    hold on
    grid on
    %plot(psal_campagne,pres_campagne,'+r')
    %plot(psal_argo_optim,pres_argo_corr,'+b')
    %plot(psal_argo_median,pres_argo_corr,'+g')
    plot(difference_psal_theta_optim,pres_argo_corr,'+r')
    plot(difference_psal_theta_median,pres_argo_corr,'+g')
    
    maxval=max(abs(difference_psal_theta_median(pres_argo>str2num(minDEPTH))));
    if isempty(maxval)||isnan(maxval)
         maxval=max(abs(difference_psal_theta_median(pres_argo>str2num('1000'))));
     end
	if isempty(maxval)||isnan(maxval)
         maxval=max(abs(difference_psal_theta_median(pres_argo>str2num('500'))));
     end
    %set(gca,'XLim',[-maxval-0.001,maxval+0.001])
    set(gca,'XLim',[-0.01,0.01])
    set(gca,'Ydir','reverse')
    grid on
    box on
    title({['Salinity . Float ' floatnum ', Cycle ' nocyclstr ]},'Fontsize',14)
    xlabel('Salinity ','Fontsize',12)
    legend({['Modified profile: optimal CPCOR (', num2str(cnew(1)), ') , optimal offset (', num2str(theoffsset_op), ')'],['Modified profile: median CPCOR (', num2str(CPcor_fixed), ') , optimal offset (', num2str(theoffsset_med), ')']},'location','SouthOutside','Fontsize',14,'FontWeight','bold')
    %legend({'Reference profile from the shipboard CTD', ['Modified profile: optimal CPCOR (', num2str(cnew(1)), ')' ],['Modified profile: median CPCOR (', num2str(CPcor_fixed), ')']},'location','SouthOutside','Fontsize',14,'FontWeight','bold')

    ylabel('pressure','Fontsize',12)
    
 
    
    set(gcf,'papertype','A4','paperunits','centimeters','paperorientation','portrait','paperposition',[.25,.75,20,25]);
    eval(['print -dpdf ' dir_enregistre '/CPCOR_analysis_median_' floatnum '_' nocyclstr   '_new.pdf']);
   % eval(['print -dpdf ' 'test.pdf']);
    fprintf(fw,'%s\n',[floatnum ', ' nocyclstr ', ' minDEPTH  ', ' num2str(cnew(1)) ', ' num2str(theoffsset_nom) ', ' num2str(theoffsset_op) ', ' num2str(theoffsset_med) ', ' num2str(theoffsset_pres) ', ' num2str(cnew(2)) ])
    %fprintf(fw,'%s\n',[floatnum ', ' nocyclstr   ', ' num2str(cnew(1)) ', ' num2str(theoffsset_op)])
    end
    %close all
 end
 keyboard
 theres=NaN*ones(117,131);
 test_cp=[-3e-7:3e-9:0.5e-7];
 test_M=[0.9993:1e-5:1.0006];
 for ir=[1:117]
 ir
   for jr=[1:131]
       theres(ir,jr)=myfun([test_cp(ir),test_M(jr)]);
    end
end
figure
pcolor(test_M,test_cp,theres)
colorbar
print -dpdf residus.pdf
theres2=NaN*ones(117,131);
%test_cp=[cnew(1)-:3e-9:0.5e-7];
test_M=[0.9998:3e-6:1.0002];
 for ir=[1:117]
 ir
   for jr=[1:134]
       theres(ir,jr)=myfun([test_cp(ir),test_M(jr)]);
    end
end
for jr=[1:131]
     [themin(jr),indic(jr)]=min(theres(:,jr));
end
for ir=[1:117]
     [themin(ir),indic(ir)]=min(theres(ir,:));
end
figure
plot(test_M,themin)
colorbar
print -dpdf residus_fond_valle.pdf
 	 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = myfun(inp)
global CPcor
global CTcor
global temp_argo pres_argo psal_argo lon_argo lat_argo psal_campagne ct_campagne pres_campagne
global minDEPTH
global press_offset
global CPcor_fixed
global M0 

% back out the float raw conservative conductivities using the nominal CTD
% values from SBE.

a1 = (1 + CTcor.*temp_argo + CPcor.*pres_argo);


cond_argo = gsw_C_from_SP(psal_argo,temp_argo,pres_argo);
cond_argo_raw = cond_argo.*a1;

pres_argo_corr=pres_argo + press_offset;

psal_argo_corr = gsw_SP_from_C(cond_argo,temp_argo,pres_argo_corr);

%cco_argo_raw = cco_argo.*a1;
% now calculate the optimized float conservative conductivities cco_opt
% using the (after the optimization iterations have converged)
% observationally estimated CPcor_new value and a multiplicative
% calibration value M.

if length(inp)>1
    CPcor_new=inp(1);
    M=inp(2);
else
    CPcor_new= CPcor_fixed;
    M=inp(1);
end
%b1 = (1 + CTcor.*temp_argo + (CPcor_new*Cpscale +Cp0).*pres_argo_corr);
b1 = (1 + CTcor.*temp_argo + (CPcor_new).*pres_argo_corr);
%cond_opt = (M*Mscale+M0)*cond_argo_raw./b1;
M0
cond_opt = M.*M0.*cond_argo_raw./b1;


psal_argo_opt= gsw_SP_from_C(cond_opt,temp_argo,pres_argo_corr);
sa_argo = gsw_SA_from_SP(psal_argo_opt,pres_argo_corr,lon_argo,lat_argo);
ct_argo = gsw_CT_from_t(sa_argo,temp_argo,pres_argo_corr);
cco_argo_opt = gsw_C_from_SP(psal_argo_opt,ct_argo,0*pres_argo_corr);

% interpolation on float theta's levels
[psal_i_campagne,pres_i_campagne] = interp_climatology(psal_campagne',ct_campagne',pres_campagne',psal_argo_opt,ct_argo,pres_argo_corr); % routine OW
cco_i_campagne = gsw_C_from_SP(psal_i_campagne,ct_argo',0*pres_i_campagne);
% Difference de salinite sur les niveaux theta

difference_psal_theta = [psal_argo_opt-psal_i_campagne'];
difference_pres_theta = [pres_argo_corr-pres_i_campagne'];

% compute the residual of the optimized float conservative conductivity and
% the reference station conductivity
%keyboard
%dcco_opt =log(cond_opt)-log(cco_i_campagne');
dcco_opt = cco_argo_opt-cco_i_campagne';
%  % moyenne glissante sur 1000db
%  pres_limup=[500:100:3500];
%  pres_liminf=[1500:100:4500];
%
%  for iij=1:length(pres_limup)
%  kkl=find(pres_argo>=pres_limup(iij)&pres_argo<=pres_liminf(iij));
%  if isempty(kkl)==0
%  dcco_opt_gliss(iij)=meanoutnan(dcco_opt(kkl));
%  else
%  dcco_opt_gliss(iij)=NaN;
%  end
%  end
%  dcco_opt=dcco_opt_gliss;
ii=find(isfinite(dcco_opt)==1 & pres_argo_corr > str2num(minDEPTH) & abs(difference_pres_theta)<1500);
if isempty(ii)
ii=find(isfinite(dcco_opt)==1 & pres_argo_corr > 1000 & abs(difference_pres_theta)<1500);
end
if isempty(ii)
ii=find(isfinite(dcco_opt)==1 & pres_argo_corr > 500 & abs(difference_pres_theta)<1500);
end
%ii=find(isfinite(dcco_opt)==1 &choosen_depth'==1 & abs(difference_pres_theta)<1500);
%ii=find(isfinite(dcco_opt)==1);
% sous echantillonne p pour ne pas sur-representer la surface
isok=1;
ide(1)=1;
hi=1;

while isok==1
    hi=hi+1;
    a=diff(pres_argo_corr(ii));
    b=cumsum(a(ide(hi-1):end));
    [~,h]=min(abs(b-a(end)));
    p=ide(hi-1)+h;
    if p<length(pres_argo_corr(ii))
        ide(hi)=p;
    else
        isok=0;
    end
end

% minimize of the sum of the squares of the residual for pressures
% exceeding the threshold value.  This is a bit heavily weighted toward the
% upper water column since the sampling is finer there. One could put the
% data on a pressure grid with more uniform intervals prior to
% minimization. One could also use an L1 Norm rather than an L2 norm if one
% wanted to get fancy.
restomin=dcco_opt(ii(ide));
%restomin=dcco_opt(ii);
std_res=std(restomin);mean_res=mean(restomin);
io=find(restomin<mean_res+3*std_res&restomin>mean_res-3*std_res);
% on n'essaie pas de minimiser les outliers

res=sum(restomin(io).^2); % L2 Norm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psal_argo_new = change_cpcor(cnew)

global CPcor
global CTcor
global temp_argo pres_argo cco_argo ct_argo psal_argo
global press_offset

% back out the float raw conductivities using the nominal CTD
% values from SBE.

pres_argo_corr = pres_argo + press_offset;

cond_argo = gsw_C_from_SP(psal_argo,temp_argo,pres_argo);

a1 = (1 + CTcor.*temp_argo + CPcor.*pres_argo);
cond_argo_raw = cond_argo.*a1;

% now calculate the optimized float conservative conductivities cco_opt
% using the (after the optimization iterations have converged)
% observationally estimated CPcor_new value and a multiplicative
% calibration value M.

CPcor_new=cnew(1);
%CTcor_new=inp(2);
%start here
M=cnew(2);
%M=1;
b1 = (1 + CTcor.*temp_argo + CPcor_new.*pres_argo_corr);

cond_argo_new=M*cond_argo_raw./b1;

% compute the corresponding psal
psal_argo_new=gsw_SP_from_C(cond_argo_new,temp_argo,pres_argo_corr);

%-------------------------------------------------------------------------------
%
% mltnc_argofile_create    - Ecriture des variables
%                            du fichier multistations au format NetCDF 		
%
%-------------------------------------------------------------------------------
%  version:
%  --------
%  1.01   Creation              			15/05/03  E. Autret
%       
% Modification Carole GRIT / Fevrier 2004 pour Etude SPMW
% 
%
%-------------------------------------------------------------------------------
function mltnc_argofile_wvar_pres(ncfile_name,data_type,format_version,...
    reference_date_time,project_name,...
    pi_name,platform_number,...
    cycle_number,direction,dc_reference,data_state_indicator,...
    data_mode,inst_reference,wmo_inst_type,juld,latitude,longitude,...
    pres_std,temp_std,psal_std,teta_std,sig0_std,...
    vorp_std,qc_prs,qc_tmp,qc_sal,qc_tet,qc_si0,qc_vrp,...
    pres_std_err,temp_std_err,tet_std_err,...
    psal_std_err,sigi_std_err,vorp_std_err);
    
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

nc=netcdf(ncfile_name,'write');

FILLVALUE=99999;
inok=find(~finite(pres_std));
if ~isempty(inok)
  pres_std(inok)=FILLVALUE;
  pres_std_err(inok)=FILLVALUE;
  qc_prs(inok)='0';
end

inok=find(~finite(temp_std));
if ~isempty(inok)
  temp_std(inok)=FILLVALUE;
  temp_std_err(inok)=FILLVALUE;
  qc_tmp(inok)='0';
end

inok=find(~finite(psal_std));
if ~isempty(inok)
  psal_std(inok)=FILLVALUE;
  psal_std_err(inok)=FILLVALUE;
  qc_sal(inok)='0';
end

inok=find(~finite(teta_std));
if ~isempty(inok)
  teta_std(inok)=FILLVALUE;
  tet_std_err(inok)=FILLVALUE;
  qc_tet(inok)='0';
end

inok=find(~finite(sig0_std));
if ~isempty(inok)
  sig0_std(inok)=FILLVALUE;
  sigi_std_err(inok)=FILLVALUE;
  qc_si0(inok)='0';
end

inok=find(~finite(vorp_std));
if ~isempty(inok)
  vorp_std(inok)=FILLVALUE;
  vorp_std_err(inok)=FILLVALUE;
  qc_vrp(inok)='0';
end

%===============
% Informations
%===============

nc{'DATA_TYPE'}(:) = data_type;

nc{'FORMAT_VERSION'}(:) = format_version;

nc{'REFERENCE_DATE_TIME'}(:) = reference_date_time;

%===============



nc{'PROJECT_NAME'}(:,:) = project_name;
nc{'PI_NAME'}(:,:) = pi_name;
nc{'PLATFORM_NUMBER'}(:,:) = platform_number;
nc{'CYCLE_NUMBER'}(:) = cycle_number;
nc{'DIRECTION'}(:) = direction;
nc{'DC_REFERENCE'}(:,:) = dc_reference;
nc{'DATA_STATE_INDICATOR'}(:,:) = data_state_indicator;
nc{'DATA_MODE'}(:) = data_mode;
nc{'INST_REFERENCE'}(:,:) = inst_reference;
nc{'WMO_INST_TYPE'}(:,:) = wmo_inst_type;
nc{'JULD'}(:) = juld;
nc{'LATITUDE'}(:) = latitude;
nc{'LONGITUDE'}(:) = longitude;

nc{'PRES'}(:,:) = pres_std;
nc{'PRES_QC'}(:,:) = qc_prs;
nc{'PRES_ERR'}(:,:) = pres_std_err;

nc{'TEMP'}(:,:) = temp_std;
nc{'TEMP_QC'}(:,:) = qc_tmp;
nc{'TEMP_ERR'}(:,:) = temp_std_err;

nc{'PSAL'}(:,:) = psal_std;
nc{'PSAL_QC'}(:,:) = qc_sal;
nc{'PSAL_ERR'}(:,:) = psal_std_err;

nc{'TPOT'}(:,:) = teta_std;
nc{'TPOT_QC'}(:,:) = qc_tet;
nc{'TPOT_ERR'}(:,:) = tet_std_err;

nc{'SIG0'}(:,:) = sig0_std;
nc{'SIG0_QC'}(:,:) = qc_si0;
nc{'SIG0_ERR'}(:,:) = sigi_std_err;

nc{'VORP'}(:,:) = vorp_std;
nc{'VORP_QC'}(:,:) = qc_vrp;
nc{'VORP_ERR'}(:,:) = vorp_std_err;

%nc{'VORP_LISS'}(:,:) = vorp_std_lissee;
%nc{'VORP_LISS_ERR'}(:,:) = vorp_liss_std_err;


endef(nc);
close(nc);







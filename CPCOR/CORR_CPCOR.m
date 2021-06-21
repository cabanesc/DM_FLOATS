% -========================================================
%   USAGE :   CORR_CPCOR(floatname,dacname,varargin)
%   PURPOSE : correct Cpcor value in R files
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char) e.g.  'coriolis'
%   OPTIONNAL INPUT : 
%    'NEW_CPCOR'     (float)    '13.5e-8' (default):  CPCor values used to correct float salinity data
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created C. Cabanes - 2020
% ========================================================
function CORR_CPCOR(floatname,dacname,varargin)
close all

init_path

if iscell(floatname)==0
floatname=cellstr(floatname);
end
if iscell(dacname)==0
dacname=cellstr(dacname)
end

CONFIG=load_configuration('config.txt');

n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end


f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);

PARAM.NEW_CPCOR=-13.5e-8;
if isfield(s,'NEW_CPCOR')==1;PARAM.NEW_CPCOR=s.NEW_CPCOR;end;

disp(' ')
disp(['NEW CPCOR is ' num2str(PARAM.NEW_CPCOR)])


% find Core files.
for ifloat=1:length(floatname)
	thefloatname=floatname{ifloat};
	thedacname=dacname{ifloat};
	IncludeDescProf=1;
	[file_list] = select_float_files_on_ftp(thefloatname,thedacname,CONFIG.DIR_FTP,'C',IncludeDescProf);
	% on remplit les champs adjusted (sauf pour pres s'i y a deja eu un ajustement)
	FILENAME_TECH = [CONFIG.DIR_FTP thedacname '/' thefloatname '/' thefloatname '_tech.nc'];
	T = read_netcdf_allthefile(FILENAME_TECH);
	% trouve les pressions de surface
	isurf=find(findstr_tab(cellstr(T.technical_parameter_name.data),'PRES_SurfaceOffsetCorrectedNotResetNegative_1cBarResolution_dbar')==1);
	ivolt=find(findstr_tab(cellstr(T.technical_parameter_name.data),'VOLTAGE_BatteryPumpStartProfile_volts')==1);
	ipump=find(findstr_tab(cellstr(T.technical_parameter_name.data),'NUMBER_PumpActionsDuringAscentToSurface_COUNT')==1);
    ivalve=find(findstr_tab(cellstr(T.technical_parameter_name.data),'NUMBER_ValveActionsDuringDescentToPark_COUNT')==1);
	cycl=T.cycle_number.data(isurf);
	cyc1=(cycl==1);
	
	figure
	subplot(3,1,1)
	plot(T.cycle_number.data(isurf),str2num(T.technical_parameter_value.data(isurf,:)),'+-')
	hold on
	box on
	grid on
	xlabel('cycle_number','interpreter','none')
	ylabel('dbar')
	title('PRES_SurfaceOffsetCorrectedNotResetNegative_1cBarResolution_dbar','interpreter','none')
	subplot(3,1,2)
	plot(T.cycle_number.data(ivolt),str2num(T.technical_parameter_value.data(ivolt,:)),'+-')
	hold on
	box on
	grid on
	xlabel('cycle_number','interpreter','none')
	ylabel('volt')
	title('VOLTAGE_BatteryPumpStartProfile_volts','interpreter','none')
	subplot(3,1,3)
	plot(T.cycle_number.data(ipump),str2num(T.technical_parameter_value.data(ipump,:)),'+-')
	hold on
	plot(T.cycle_number.data(ivalve),str2num(T.technical_parameter_value.data(ivalve,:)),'o-r')
	hold on
	box on
	grid on
	xlabel('cycle_number','interpreter','none')
	ylabel('COUNT')
	title({'NUMBER_PumpActionsDuringAscentToSurface_COUNT (+)'; 'NUMBER_ValveActionsDuringDescentToPark_COUNT (o)'},'interpreter','none')
% 	fileout=[CONFIG.DIR_PLOT 'preliminaire/' thefloatname '/' thefloatname '_surface_pres.pdf'];
% 	if exist([CONFIG.DIR_PLOT 'preliminaire/' thefloatname '/'])==0
% 	mkdir([CONFIG.DIR_PLOT 'preliminaire/' thefloatname '/'])
% 	end
% 	%keyboard
% 	set(gcf,'Position',[1 46 660 800],'paperPositionMode','auto');
% 	print(fileout,'-dpdf')
	
	
	
	for ifiles=1:length(file_list)

		file_name = [,CONFIG.DIR_FTP thedacname '/' thefloatname '/profiles/' file_list{ifiles} ];
		[F,Dim,G]=read_netcdf_allthefile(file_name);
		F = replace_fill_bynan(F);
		Forig=F;

		if F.cycle_number.data==1
		press_offset= str2num(T.technical_parameter_value.data(isurf(cyc1),:));
		else
		press_offset=0;
		end


		F.temp_adjusted.data=F.temp.data;
		F.temp_adjusted_qc.data=F.temp_qc.data;

		test = check_isfillval_prof(F,'pres_adjusted_error');
		if test.pres_adjusted_error==1
		F.pres_adjusted.data=F.pres.data+press_offset;
		F.pres_adjusted_qc.data=F.pres_qc.data;
		end

		F.psal_adjusted_qc.data = F.psal_qc.data;



		cnew=[PARAM.NEW_CPCOR,1];
		psal_argo_new = change_cpcor(F,cnew);
		
		meancor(ifiles)=meanoutnan(psal_argo_new(1,F.pres.data(1,:)>=500)-F.psal.data(1,F.pres.data(1,:)>=500));
		maxpres(ifiles)=max(F.pres.data(1,:));

		isprimary = findstr_tab(F.vertical_sampling_scheme.data,'Primary sampling');
		F.psal_adjusted.data(isprimary,:) = psal_argo_new(isprimary,:);

		n_calib=1;
		theparameters = strtrim(squeeze(F.parameter.data(isprimary,n_calib,:,:)));
		ind_psal = ismember(cellstr(theparameters),'PSAL');

		thestr=(['CCo_ADJUSTED=CCo*(1+delta*TEMP+CpcorSBE*PRES)/((1+delta*TEMP+Cpcor*PRES);  delta=+3.25e-6; CPcorSBE=-9.57e-08.']);
		l_thestr=length(thestr);
		F.scientific_calib_equation.data(isprimary,n_calib,ind_psal,:) = F.scientific_calib_equation.FillValue_;
		F.scientific_calib_equation.data(isprimary,n_calib,ind_psal,1:l_thestr) = thestr;
		thestr=(['CPcor= -12.5e-8']);
		l_thestr=length(thestr);
		F.scientific_calib_coefficient.data(isprimary,n_calib,ind_psal,:) = F.scientific_calib_coefficient.FillValue_;
		F.scientific_calib_coefficient.data(isprimary,n_calib,ind_psal,1:l_thestr) = thestr;   
		thestr=(['Recomputation of the conservative conductivity CCo using a different CPcor value from that provided by Seabird']);
		l_thestr=length(thestr);
		F.scientific_calib_comment.data(isprimary,n_calib,ind_psal,:) = F.scientific_calib_comment.FillValue_;
		F.scientific_calib_comment.data(isprimary,n_calib,ind_psal,1:l_thestr) = thestr;

		if F.cycle_number.data==1
				ind_pres = ismember(cellstr(theparameters),'PRES');
				
				thestr=(['PRES_ADJUSTED(cycle 1)=PRES (cycle 1)-Surface Pressure(cycle 1)']);
				l_thestr=length(thestr);
				F.scientific_calib_equation.data(isprimary,n_calib,ind_psal,:) = F.scientific_calib_equation.FillValue_;
				F.scientific_calib_equation.data(isprimary,n_calib,ind_psal,1:l_thestr) = thestr;
				thestr=(['Surface Pressure =' num2str(-press_offset)]);
				l_thestr=length(thestr);
				F.scientific_calib_coefficient.data(isprimary,n_calib,ind_psal,:) = F.scientific_calib_coefficient.FillValue_;
				F.scientific_calib_coefficient.data(isprimary,n_calib,ind_psal,1:l_thestr) = thestr;   
				thestr=(['This float is autocorrecting but press_offset is  0 at cycle 1. Instead ,we used press_offset =-Surface Pressure(cycle 1)']);
				l_thestr=length(thestr);
				F.scientific_calib_comment.data(isprimary,n_calib,ind_psal,:) = F.scientific_calib_comment.FillValue_;
				F.scientific_calib_comment.data(isprimary,n_calib,ind_psal,1:l_thestr) = thestr;
		end


		Fres=F;
		F=replace_nan_byfill(F);
		%create_netcdf_allthefile(F,Dim,'test.nc',G)
		create_netcdf_allthefile(F,Dim,file_name,G)


	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psal_argo_new = change_cpcor(F,cnew)

% nominal CT values from SBE.
%--------------------------------------------------------------------------
CPcor = -9.5700E-8;
CTcor =  3.2500E-6;

global temp_argo pres_argo cco_argo ct_argo psal_argo
global press_offset

% back out the float raw conductivities using the nominal CTD
% values from SBE.
pres_argo = F.pres.data;
pres_argo_corr = F.pres_adjusted.data;
psal_argo=F.psal.data;
temp_argo=F.temp.data;

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


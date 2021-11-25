% -========================================================
%   USAGE :   CORR_CPCOR(floatname,dacname,varargin)
%   PURPOSE : compare first descending or ascending profile to ctd reference made at float launch
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char) e.g.  'coriolis'
%   OPTIONNAL INPUT : 
%    'NEW_CPCOR'     (float)    '-13.5e-8' (default): median CPCor values used to correct float salinity data
%    'NEW_M'         (float)    '1' (default) : multiplicative factor for
%    conservative conductivity
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created C. Cabanes - 2020
% ========================================================
function CORR_CPCOR(floatname,dacname,config_campaign_file,varargin)
%close all

%init_path

if iscell(floatname)==0
floatname=cellstr(floatname);
end
if iscell(dacname)==0
dacname=cellstr(dacname)
end

CONFIG=load_configuration('config.txt');
CAMPAIGN=load_configuration(config_campaign_file);

n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end


f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);

CPcor_SBE = -9.5700E-8;
CTcor =  3.2500E-6;

PARAM.NEW_CPCOR=-13.5e-8;
PARAM.NEW_M=1;
if isfield(s,'NEW_CPCOR')==1;PARAM.NEW_CPCOR=s.NEW_CPCOR;end;
if isfield(s,'NEW_M')==1;PARAM.NEW_M=s.NEW_M;end;
disp(' ')
disp(['NEW CPCOR is ' num2str(PARAM.NEW_CPCOR)])
disp(['NEW M value is ' num2str(PARAM.NEW_M)])

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
	
% 	figure
% 	subplot(3,1,1)
% 	plot(T.cycle_number.data(isurf),str2num(T.technical_parameter_value.data(isurf,:)),'+-')
% 	hold on
% 	box on
% 	grid on
% 	xlabel('cycle_number','interpreter','none')
% 	ylabel('dbar')
% 	title('PRES_SurfaceOffsetCorrectedNotResetNegative_1cBarResolution_dbar','interpreter','none')
% 	subplot(3,1,2)
% 	plot(T.cycle_number.data(ivolt),str2num(T.technical_parameter_value.data(ivolt,:)),'+-')
% 	hold on
% 	box on
% 	grid on
% 	xlabel('cycle_number','interpreter','none')
% 	ylabel('volt')
% 	title('VOLTAGE_BatteryPumpStartProfile_volts','interpreter','none')
% 	subplot(3,1,3)
% 	plot(T.cycle_number.data(ipump),str2num(T.technical_parameter_value.data(ipump,:)),'+-')
% 	hold on
% 	plot(T.cycle_number.data(ivalve),str2num(T.technical_parameter_value.data(ivalve,:)),'o-r')
% 	hold on
% 	box on
% 	grid on
% 	xlabel('cycle_number','interpreter','none')
% 	ylabel('COUNT')
% 	title({'NUMBER_PumpActionsDuringAscentToSurface_COUNT (+)'; 'NUMBER_ValveActionsDuringDescentToPark_COUNT (o)'},'interpreter','none')
% 	fileout=[CONFIG.DIR_PLOT 'preliminaire/' thefloatname '/' thefloatname '_surface_pres.pdf'];
% 	if exist([CONFIG.DIR_PLOT 'preliminaire/' thefloatname '/'])==0
% 	mkdir([CONFIG.DIR_PLOT 'preliminaire/' thefloatname '/'])
% 	end
% 	%keyboard
% 	set(gcf,'Position',[1 46 660 800],'paperPositionMode','auto');
% 	print(fileout,'-dpdf')
	fileout=[CONFIG.DIR_PLOT 'verif_profil1/' thefloatname '/cpcor_value_applied.txt'];
	if exist([CONFIG.DIR_PLOT 'verif_profil1/' thefloatname '/'])==0
	mkdir([CONFIG.DIR_PLOT 'verif_profil1/' thefloatname '/'])
	end
	fw1=fopen(fileout,'w');
	fprintf(fw1,'%s\n', num2str(PARAM.NEW_CPCOR));
	fclose(fw1)
	
	
	for ifiles=1:length(file_list)

		file_name = [,CONFIG.DIR_FTP thedacname '/' thefloatname '/profiles/' file_list{ifiles} ];
		[F,Dim,G]=read_netcdf_allthefile(file_name);
		F = replace_fill_bynan(F);
		Forig=F;

		if F.cycle_number.data==1
		%press_offset= str2num(T.technical_parameter_value.data(isurf(cyc1),:)); 
		% Update 04/2021: ne tient pas compte de l'offset normalement c'est corrigé (voir le manuel utilisateur des deep aRG0
		% Le pres offset est mesuré avant la descente, est gardé en memoire et utilisé pour corriger les pressions au cours du cycle. transmis 10 jours plus tard.
		% Si correction de pression cycle 1A, on remarque que la plupart des pressions near surface sont négatives!!! pas cohérent.
		press_offset=0;
		else
		press_offset=0;
		end
        
		isprimary = find(findstr_tab(F.vertical_sampling_scheme.data,'Primary sampling'));
		notisprimary = find(findstr_tab(F.vertical_sampling_scheme.data,'Near-surface sampling'));
        F.data_mode.data(notisprimary)='A';
        %keyboard
		F.temp_adjusted.data=F.temp.data;
		F.temp_adjusted_qc.data=F.temp_qc.data;
		
		
	    F.temp_adjusted_error.data=0.002*ones(Dim.n_prof.dimlength,Dim.n_levels.dimlength);
	    F.temp_adjusted_error.data(isnan(F.temp_adjusted.data))=NaN;
	

		test = check_isfillval_prof(F,'pres_adjusted_error');
		%if test.pres_adjusted_error==1 % pour les arvor/provor pas de correction de pression RT cycle 1
		F.pres_adjusted.data=F.pres.data+press_offset;
		F.pres_adjusted_qc.data=F.pres_qc.data;
		%end
        F.pres_adjusted_error.data=(2.5/6000)*F.pres_adjusted.data+2; % improved error for deep floats.
		
		F.psal_adjusted_qc.data = F.psal_qc.data;
        F.psal_adjusted_error.data=0.004*ones(Dim.n_prof.dimlength,Dim.n_levels.dimlength);
	    F.psal_adjusted_error.data(isnan(F.temp_adjusted.data))=NaN;

		cnew=[PARAM.NEW_CPCOR,PARAM.NEW_M];
		psal_argo_new = change_cpcor(F,cnew);

		F.psal_adjusted.data = psal_argo_new;

		n_calib=1;
		theparameters = strtrim(squeeze(F.parameter.data(isprimary,n_calib,:,:)));
		ind_psal = ismember(cellstr(theparameters),'PSAL');
        
		
		thestr=(['new conductivity = original conductivity * (1 + delta*TEMP + CPcor_SBE*PRES) / (1 + delta*TEMP_ADJUSTED + CPcor_new*PRES_ADJUSTED)']);
		l_thestr=length(thestr);
		F.scientific_calib_equation.data(isprimary,n_calib,ind_psal,:) = F.scientific_calib_equation.FillValue_;
		F.scientific_calib_equation.data(isprimary,n_calib,ind_psal,1:l_thestr) = thestr;
		if ~isempty (notisprimary)
		F.scientific_calib_equation.data(notisprimary,n_calib,ind_psal,:) = F.scientific_calib_equation.FillValue_;
		F.scientific_calib_equation.data(notisprimary,n_calib,ind_psal,1:l_thestr) = thestr;
		end
		
		
        thestr=(['CPcor_new = ' num2str(PARAM.NEW_CPCOR) '; CPcor_SBE = ' num2str(CPcor_SBE) '; delta = ' num2str(CTcor)]);		l_thestr=length(thestr);
		F.scientific_calib_coefficient.data(isprimary,n_calib,ind_psal,:) = F.scientific_calib_coefficient.FillValue_;
		F.scientific_calib_coefficient.data(isprimary,n_calib,ind_psal,1:l_thestr) = thestr;   
		if ~isempty (notisprimary)
		F.scientific_calib_coefficient.data(notisprimary,n_calib,ind_psal,:) = F.scientific_calib_coefficient.FillValue_;
		F.scientific_calib_coefficient.data(notisprimary,n_calib,ind_psal,1:l_thestr) = thestr; 
		end
		
		thestr=(['New conductivity computed by using a different CPcor value from that provided by Sea-Bird.']);
		l_thestr=length(thestr);
		F.scientific_calib_comment.data(isprimary,n_calib,ind_psal,:) = F.scientific_calib_comment.FillValue_;
		F.scientific_calib_comment.data(isprimary,n_calib,ind_psal,1:l_thestr) = thestr;

		thedate=datestr(now,'yyyymmddHHMMSS');
		F.scientific_calib_date.data(isprimary,n_calib,ind_psal,:)=thedate;

		if ~isempty (notisprimary)
		F.scientific_calib_comment.data(notisprimary,n_calib,ind_psal,:) = F.scientific_calib_comment.FillValue_;
		F.scientific_calib_comment.data(notisprimary,n_calib,ind_psal,1:l_thestr) = thestr;
		F.scientific_calib_date.data(notisprimary,n_calib,ind_psal,:)=thedate;
		end
       

		ind_pres = ismember(cellstr(theparameters),'PRES');

		% if F.cycle_number.data==1
				
				% thestr=(['PRES_ADJUSTED(cycle 1)=PRES (cycle 1)- Surface Pressure(cycle 1)']);
				% l_thestr=length(thestr);
				% F.scientific_calib_equation.data(isprimary,n_calib,ind_pres,:) = F.scientific_calib_equation.FillValue_;
				% F.scientific_calib_equation.data(isprimary,n_calib,ind_pres,1:l_thestr) = thestr;
				% if ~isempty (notisprimary)
				% F.scientific_calib_equation.data(notisprimary,n_calib,ind_pres,:) = F.scientific_calib_equation.FillValue_;
				% F.scientific_calib_equation.data(notisprimary,n_calib,ind_pres,1:l_thestr) = thestr;
				% end
				% thestr=(['Surface Pressure =' num2str(-press_offset)]);
				% l_thestr=length(thestr);
				% F.scientific_calib_coefficient.data(isprimary,n_calib,ind_pres,:) = F.scientific_calib_coefficient.FillValue_;
				% F.scientific_calib_coefficient.data(isprimary,n_calib,ind_pres,1:l_thestr) = thestr;  
				% if ~isempty (notisprimary)
				% F.scientific_calib_coefficient.data(notisprimary,n_calib,ind_pres,:) = F.scientific_calib_coefficient.FillValue_;
				% F.scientific_calib_coefficient.data(notisprimary,n_calib,ind_pres,1:l_thestr) = thestr;  
				% end
				% thestr=(['This float is autocorrecting but not at cycle 1. We corrected the pressure with Surface Pressure(cycle 1)']);
				% l_thestr=length(thestr);
				% F.scientific_calib_comment.data(isprimary,n_calib,ind_pres,:) = F.scientific_calib_comment.FillValue_;
				% F.scientific_calib_comment.data(isprimary,n_calib,ind_pres,1:l_thestr) = thestr;
				% F.scientific_calib_date.data(isprimary,n_calib,ind_pres,:)=thedate;
				% if ~isempty (notisprimary)
				% F.scientific_calib_comment.data(notisprimary,n_calib,ind_pres,:) = F.scientific_calib_comment.FillValue_;
				% F.scientific_calib_comment.data(notisprimary,n_calib,ind_pres,1:l_thestr) = thestr;
				% F.scientific_calib_date.data(notisprimary,n_calib,ind_pres,:)=thedate;
				% end
		% else
		
		        thestr=(['PRES_ADJUSTED = PRES']);
				l_thestr=length(thestr);
				F.scientific_calib_equation.data(isprimary,n_calib,ind_pres,:) = F.scientific_calib_equation.FillValue_;
				F.scientific_calib_equation.data(isprimary,n_calib,ind_pres,1:l_thestr) = thestr;
				if ~isempty (notisprimary)
				F.scientific_calib_equation.data(notisprimary,n_calib,ind_pres,:) = F.scientific_calib_equation.FillValue_;
				F.scientific_calib_equation.data(notisprimary,n_calib,ind_pres,1:l_thestr) = thestr;
				end
                thestr=(['none']);
				l_thestr=length(thestr);
				F.scientific_calib_coefficient.data(isprimary,n_calib,ind_pres,:) = F.scientific_calib_coefficient.FillValue_;
				F.scientific_calib_coefficient.data(isprimary,n_calib,ind_pres,1:l_thestr) = thestr; 
				if ~isempty (notisprimary)
				F.scientific_calib_coefficient.data(notisprimary,n_calib,ind_pres,:) = F.scientific_calib_coefficient.FillValue_;
				F.scientific_calib_coefficient.data(notisprimary,n_calib,ind_pres,1:l_thestr) = thestr; 
				end
				thestr=(['This float is autocorrecting pressure. Data is good within the specified error.']);
				l_thestr=length(thestr);
				F.scientific_calib_comment.data(isprimary,n_calib,ind_pres,:) = F.scientific_calib_comment.FillValue_;
				F.scientific_calib_comment.data(isprimary,n_calib,ind_pres,1:l_thestr) = thestr;
				F.scientific_calib_date.data(isprimary,n_calib,ind_pres,:)=thedate;
				if ~isempty (notisprimary)
				F.scientific_calib_comment.data(notisprimary,n_calib,ind_pres,:) = F.scientific_calib_comment.FillValue_;
				F.scientific_calib_comment.data(notisprimary,n_calib,ind_pres,1:l_thestr) = thestr;
				F.scientific_calib_date.data(notisprimary,n_calib,ind_pres,:)=thedate;
				end
		% end
		
		
		ind_temp = ismember(cellstr(theparameters),'TEMP');
        thestr=(['TEMP_ADJUSTED = TEMP']);
		l_thestr=length(thestr);
		F.scientific_calib_equation.data(isprimary,n_calib,ind_temp,:) = F.scientific_calib_equation.FillValue_;
		F.scientific_calib_equation.data(isprimary,n_calib,ind_temp,1:l_thestr) = thestr;
		if ~isempty (notisprimary)
        F.scientific_calib_equation.data(notisprimary,n_calib,ind_temp,:) = F.scientific_calib_equation.FillValue_;
		F.scientific_calib_equation.data(notisprimary,n_calib,ind_temp,1:l_thestr) = thestr;
		end
		thestr=(['none']);
		l_thestr=length(thestr);
		F.scientific_calib_coefficient.data(isprimary,n_calib,ind_temp,:) = F.scientific_calib_coefficient.FillValue_;
		F.scientific_calib_coefficient.data(isprimary,n_calib,ind_temp,1:l_thestr) = thestr; 
		if ~isempty (notisprimary)
		F.scientific_calib_coefficient.data(notisprimary,n_calib,ind_temp,:) = F.scientific_calib_coefficient.FillValue_;
		F.scientific_calib_coefficient.data(notisprimary,n_calib,ind_temp,1:l_thestr) = thestr; 
		end
		thestr=(['Data is good within the specified error.']);
		l_thestr=length(thestr);
		F.scientific_calib_comment.data(isprimary,n_calib,ind_temp,:) = F.scientific_calib_comment.FillValue_;
		F.scientific_calib_comment.data(isprimary,n_calib,ind_temp,1:l_thestr) = thestr;
		F.scientific_calib_date.data(isprimary,n_calib,ind_temp,:)=thedate;
		if ~isempty (notisprimary)
		F.scientific_calib_comment.data(notisprimary,n_calib,ind_temp,:) = F.scientific_calib_comment.FillValue_;
		F.scientific_calib_comment.data(notisprimary,n_calib,ind_temp,1:l_thestr) = thestr;
		F.scientific_calib_date.data(notisprimary,n_calib,ind_temp,:)=thedate;
		end
		
        %================== HISTORY fields
	   %define new_hist=N_HISTORY+1
	    if PARAM.NEW_CPCOR~=-13.5e-8  % add a new history ste
 
			new_hist = Dim.n_history.dimlength+1;
			allfields = fieldnames(F);
			ii = strfind(allfields,'history_');
			is_history = find(~cellfun('isempty',ii));
			n_prof=isprimary;		
			F = check_FirstDimArray_is(F,'N_HISTORY');
			
			% initialize a new history        
			if Dim.n_history.dimlength~=0
				[F_ex,Dim_ex]=extract_profile_dim(F,Dim,'N_HISTORY',1);
				for ik = is_history'
					oneChamp =allfields{ik};
					ii=F_ex.(oneChamp).data~=F_ex.(oneChamp).FillValue_;
					F_ex.(oneChamp).data(ii)=F_ex.(oneChamp).FillValue_;
				end
				[F,Dim] = cat_profile_dim(F,F_ex,Dim,Dim_ex,'N_HISTORY');
			else
				for ik = is_history'
					oneChamp =allfields{ik};
					siz(1)=1;
					for tk=2:length(F.(oneChamp).dim)
						siz(tk) = Dim.(lower(F.(oneChamp).dim{tk})).dimlength;
					end
					F.(oneChamp).data = repmat(FLD_ex.(oneChamp).FillValue_,siz);
					Dim.n_history.dimlength=1;
				end
			end
		
			% fill HISTORY section
			institution=CONFIG.OPERATOR_INSTITUTION;
			l_in=length(institution);
			F.history_institution.data(new_hist,n_prof,1:l_in)=institution;
			
			step='ARSQ';
			F.history_step.data(new_hist,n_prof,:)=step;
			
			soft='DMCP';
			l_so=length(soft);
			F.history_software.data(new_hist,n_prof,1:l_so)=soft;
			
			
			soft_release='1.0';
			l_so_r=length(soft_release);
			F.history_software_release.data(new_hist,n_prof,1:l_so_r)=soft_release;		
			reference=CAMPAIGN.STATION;
			l_ref=length(reference);
			F.history_reference.data(new_hist,n_prof,1:l_ref)=reference;
			
			action='IP';
			l_ac=length(action);
			F.history_action.data(new_hist,n_prof,1:l_ac)=action;
			
			F.history_date.data(new_hist,n_prof,:)=thedate;	
			
			parameter='PSAL';
			l_pa=length(parameter);
			F.history_parameter.data(new_hist,n_prof,1:l_pa)=parameter;
		
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


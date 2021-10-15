% -========================================================
%   USAGE : plot_model_fit(tabfloat,tabdac,configow,POSTFIX_D,POSTFIX_S)
%   PURPOSE : superpose results of owc (calibration figure) with different
%   set_calseries parameters (usefull for deep floats)
% -----------------------------------
%   INPUT :
%    tabfloat  (char or cell of chars -size n_floatsx1)    e.g. '6900258' or {'6900258', '3901954'}
%    tabdac    (char or cell of chars -size n_floatsx1)    e.g. 'coriolis' or {'coriolis', 'bodc'}
%    configow  (floats or cell of floats -size n_floatsx1) e.g.  149       or {149,149}      % config number ow
%
%    'POSTFIX_D'      (char)      deepest zone e.g  '_1500_2000' document choices made in set_calseries.m for use_pres_gt use_pres_lt
%    'POSTFIX_S'      (char)     shallowest zone e.g '_500_1000' document choices made in set_calseries.m for use_pres_gt use_pres_lt
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created (2016) ccabanes
%
%   CALLED SUBROUTINES: 
% -------------------------------------
% 
% ========================================================
function plot_model_fit(tabfloat,tabdac,configow,POSTFIX_D,POSTFIX_S)

%init_path

close all

C = load_configuration('config.txt');
VERSION_OW = C.VERSION_OWC;
addpath([C.DIR_OWC VERSION_OW '/matlab_codes/']) % codes OW originaux
addpath('util/')
addpath('ow_codes_modif/')  % modif codes OW
addpath('ow_codes_modif/Speed/')  % modif codes OW
flt_name=tabfloat;


for THEPLOT=[1,2];

	pn_float_dir='';
	pn_float_name=num2str(flt_name);
	flt_name = deblank(num2str(flt_name));
	po_system_configuration.FLOAT_SOURCE_DIRECTORY= [C.DIR_DATA 'float_source/'];
	po_system_configuration.FLOAT_SOURCE_POSTFIX='.mat';
	po_system_configuration.FLOAT_MAPPED_DIRECTORY= [C.DIR_DATA 'float_mapped/CONFIG' num2str(configow) '/'];
	po_system_configuration.FLOAT_MAPPED_PREFIX='map_';
	po_system_configuration.FLOAT_MAPPED_POSTFIX='.mat';
	po_system_configuration.FLOAT_CALIB_DIRECTORY= [C.DIR_DATA 'float_calib/CONFIG' num2str(configow) '/' ];
	po_system_configuration.FLOAT_CALIB_PREFIX='cal_';
	po_system_configuration.FLOAT_CALSERIES_PREFIX='calseries_';

	if THEPLOT==1
	po_system_configuration.FLOAT_CALIB_POSTFIX=[POSTFIX_D '.mat'];
	po_system_configuration.FLOAT_CALSERIES_POSTFIX=[POSTFIX_D '.mat'];
	elseif THEPLOT==2
	po_system_configuration.FLOAT_CALIB_POSTFIX=[POSTFIX_S '.mat'];
	po_system_configuration.FLOAT_CALSERIES_POSTFIX=[POSTFIX_S '.mat'];
	end



	% modify pn_float_name for title if the name contains '_' ------------

	ii=find(pn_float_name=='_');
	if(isempty(ii)==0)
		title_floatname = strcat( pn_float_name(1:ii-1), '\_', pn_float_name(ii+1:length(pn_float_name)));
	else
		title_floatname = pn_float_name;
	end


	% load data from /float_source, /float_mapped, /float_calib, and config --------------

	lo_float_source_data = load(fullfile( po_system_configuration.FLOAT_SOURCE_DIRECTORY, pn_float_dir, ...
		strcat( pn_float_name, po_system_configuration.FLOAT_SOURCE_POSTFIX ) ) );

	Co.longitude.data=lo_float_source_data.LONG;
	Co.latitude.data=lo_float_source_data.LAT;


	Co=shiftEW(Co,'longitude','grwch');

	PROFILE_NO = lo_float_source_data.PROFILE_NO;
	LAT  = Co.latitude.data;
	LONG = Co.longitude.data;
	PRES = lo_float_source_data.PRES;
	TEMP = lo_float_source_data.TEMP;
	PTMP = lo_float_source_data.PTMP;
	SAL  = lo_float_source_data.SAL;

	if(isempty(find(isnan(PRES)==0))==0) % if no data exists, terminate here, no plots will be produced
		
		lo_float_mapped_data = load( fullfile( po_system_configuration.FLOAT_MAPPED_DIRECTORY, pn_float_dir, ...
			strcat( po_system_configuration.FLOAT_MAPPED_PREFIX, pn_float_name, po_system_configuration.FLOAT_MAPPED_POSTFIX ) ) ) ;
		
		mapped_sal  = lo_float_mapped_data.la_mapped_sal;
		mapsalerrors = lo_float_mapped_data.la_mapsalerrors;
		la_ptmp = lo_float_mapped_data.la_ptmp;
		selected_hist =lo_float_mapped_data.selected_hist;
		
		if(isempty(find(isnan(mapped_sal)==0))==0) % if mapping exists, terminate here, no plots will be produced
			
			lo_float_calib_data = load( fullfile( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, ...
				strcat( po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ) );
			
			cal_SAL = lo_float_calib_data.cal_SAL;
			cal_SAL_err = lo_float_calib_data.cal_SAL_err;
			pcond_factor = lo_float_calib_data.pcond_factor;
			pcond_factor_err = lo_float_calib_data.pcond_factor_err;
			
			lo_float_calseries = load( fullfile( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, ...
				strcat( po_system_configuration.FLOAT_CALSERIES_PREFIX , pn_float_name, po_system_configuration.FLOAT_CALSERIES_POSTFIX ) ) );
			
			use_theta_gt = lo_float_calseries.use_theta_gt;
			use_theta_lt = lo_float_calseries.use_theta_lt;
			use_pres_gt = lo_float_calseries.use_pres_gt;
			use_pres_lt = lo_float_calseries.use_pres_lt;
			use_percent_gt = lo_float_calseries.use_percent_gt;
			
			% load the station by station fits
			load(fullfile( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir,...
				strcat( po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ),'-regexp','^sta')
			
			
			% plot the float locations (figure 1) -----------------------
			
			%load( fullfile( po_system_configuration.CONFIG_DIRECTORY, po_system_configuration.CONFIG_COASTLINES ), 'coastdata_x', 'coastdata_y' );
			
			
			%keyboard
			[m,n] = size(PRES);
			
			
			
			
			
			
			% calibration curve (figure 3) --------------------------
			
			if(isempty(find(isnan(cal_SAL)==0))==0) % if no cal exists, terminate here, no plot will be produced
				
				Soffset=cal_SAL-SAL;
				avg_Soffset=NaN.*ones(1,n);
				avg_Soffset_err=NaN.*ones(1,n);
				Staoffset=sta_SAL-SAL;
				avg_Staoffset=NaN.*ones(1,n);
				avg_Staoffset_err=NaN.*ones(1,n);
				avg_mapped_err=NaN.*ones(1,n);
				
				[tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas( SAL, PTMP, PRES, la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);
				
				for i=1:n
					ii=[];
					ii=find(isnan(Soffset(:,i))==0);
					if ~isempty(ii)
						avg_Soffset(i)=mean(Soffset(ii,i));
						avg_Soffset_err(i)=mean(cal_SAL_err(ii,i));
					else
						avg_Soffset(i) = NaN;
						avg_Soffset_err(i) = NaN;
					end
					ii=find(isnan(Staoffset(:,i))==0);
					if ~isempty(ii)
						avg_Staoffset(i)=mean(Staoffset(ii,i));
						avg_Staoffset_err(i)=mean(sta_SAL_err(ii,i));
					else
						avg_Staoffset(i) = NaN;
						avg_Staoffset_err(i) = NaN;
					end
					ii=find(isnan(Soffset(:,i))==0);
					if ~isempty(ii)
						
						avg_mapped_err(i)=meanoutnan(mapsalerrors(ii,i));
					else
						
						avg_mapped_err(i) = NaN;
					end
				end
				
				% plot salinity time series on theta levels with the smallest S variance (figure 6) ------------
				
				unique_cal = unique(lo_float_calseries.calseries);
				n_seq = length(unique_cal);
				
				for iy=1:n_seq
					if unique_cal(iy)>0
						calindex = find(lo_float_calseries.calseries==unique_cal(iy));
						k = length(calindex);
						
						% choose 10 float theta levels to use in the piecewise linear fit --------
						
						unique_SAL = SAL(:, calindex);
						unique_PTMP = PTMP(:, calindex);
						unique_PRES = PRES(:, calindex);
						unique_la_ptmp = la_ptmp(:, calindex);
						unique_mapped_sal = mapped_sal(:, calindex);
						unique_mapsalerrors = mapsalerrors(:, calindex);
						
						ten_SAL = NaN.*ones(10,k);
						ten_PTMP = NaN.*ones(10,k);
						ten_PRES = NaN.*ones(10,k);
						ten_mapped_sal = NaN.*ones(10,k);
						ten_mapsalerrors = NaN.*ones(10,k);
						
						
						[tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas( unique_SAL, unique_PTMP, unique_PRES, unique_la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);
						%[tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas( unique_mapped_sal, unique_PTMP, unique_PRES, unique_la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);
						
						tplot=[1:length(tlevels)];
						p=length(tplot);
						Sint=NaN.*ones(p,n);
						Smap=NaN.*ones(p,n);
						Smaperr=NaN.*ones(p,n);
						Scal=NaN.*ones(p,n);
						Scalerr=NaN.*ones(p,n);
						Thetalevel_indexes=NaN.*ones(p,n);
						
						trimPRES=PRES;  % use only manually specified THETA & PRES range ---
						trimSAL=SAL;
						trimPTMP=PTMP;
						trim_mapped_sal=mapped_sal;
						trim_mapsalerrors=mapsalerrors;
						trim_cal_SAL=cal_SAL;
						trim_cal_SAL_err=cal_SAL_err;
						
						jj=find(isnan(la_ptmp)==1);
						trimPRES(jj)=NaN;
						trimSAL(jj)=NaN;
						trimPTMP(jj)=NaN;
						trim_mapped_sal(jj)=NaN;
						trim_mapsalerrors(jj)=NaN;
						trim_cal_SAL(jj)=NaN;
						trim_cal_SAL_err(jj)=NaN;
						
						if( isempty(use_theta_lt)==0 & isempty(use_theta_gt)==1 )
							jj=find(trimPTMP>use_theta_lt);
							trimPRES(jj)=NaN;
							trimSAL(jj)=NaN;
							trimPTMP(jj)=NaN;
							trim_mapped_sal(jj)=NaN;
							trim_mapsalerrors(jj)=NaN;
							trim_cal_SAL(jj)=NaN;
							trim_cal_SAL_err(jj)=NaN;
						end
						
						if( isempty(use_theta_gt)==0 & isempty(use_theta_lt)==1 )
							jj=find(trimPTMP<use_theta_gt);
							trimPRES(jj)=NaN;
							trimSAL(jj)=NaN;
							trimPTMP(jj)=NaN;
							trim_mapped_sal(jj)=NaN;
							trim_mapsalerrors(jj)=NaN;
							trim_cal_SAL(jj)=NaN;
							trim_cal_SAL_err(jj)=NaN;
						end
						
						if( isempty(use_theta_gt)==0 & isempty(use_theta_lt)==0 )
							if(use_theta_gt>use_theta_lt) %the middle band is excluded
								jj=find(trimPTMP<use_theta_gt&trimPTMP>use_theta_lt);
							else
								jj=find(trimPTMP<use_theta_gt|trimPTMP>use_theta_lt);
							end
							trimPRES(jj)=NaN;
							trimSAL(jj)=NaN;
							trimPTMP(jj)=NaN;
							trim_mapped_sal(jj)=NaN;
							trim_mapsalerrors(jj)=NaN;
							trim_cal_SAL(jj)=NaN;
							trim_cal_SAL_err(jj)=NaN;
						end
						
						if( isempty(use_pres_lt)==0 & isempty(use_pres_gt)==1 )
							jj=find(trimPRES>use_pres_lt);
							trimPRES(jj)=NaN;
							trimSAL(jj)=NaN;
							trimPTMP(jj)=NaN;
							trim_mapped_sal(jj)=NaN;
							trim_mapsalerrors(jj)=NaN;
							trim_cal_SAL(jj)=NaN;
							trim_cal_SAL_err(jj)=NaN;
						end
						
						if( isempty(use_pres_gt)==0 & isempty(use_pres_lt)==1 )
							jj=find(trimPRES<use_pres_gt);
							trimPRES(jj)=NaN;
							trimSAL(jj)=NaN;
							trimPTMP(jj)=NaN;
							trim_mapped_sal(jj)=NaN;
							trim_mapsalerrors(jj)=NaN;
							trim_cal_SAL(jj)=NaN;
							trim_cal_SAL_err(jj)=NaN;
						end
						
						if( isempty(use_pres_gt)==0 & isempty(use_pres_lt)==0 )
							if(use_pres_gt>use_pres_lt) %he middle band is excluded
								jj=find(trimPRES<use_pres_gt&trimPRES>use_pres_lt);
							else
								jj=find(trimPRES<use_pres_gt|trimPRES>use_pres_lt);
							end
							trimPRES(jj)=NaN;
							trimSAL(jj)=NaN;
							trimPTMP(jj)=NaN;
							trim_mapped_sal(jj)=NaN;
							trim_mapsalerrors(jj)=NaN;
							trim_cal_SAL(jj)=NaN;
							trim_cal_SAL_err(jj)=NaN;
						end
						
						for i=1:n
							for j=tplot
								if(tlevels(j)<max(trimPTMP(:,i))&tlevels(j)>min(trimPTMP(:,i)))
									diffTheta = abs(trimPTMP(:,i)-tlevels(j));
									if isempty(find(~isnan(diffTheta)))
										Thetalevel_indexes(j,i) = NaN;
									else
										Thetalevel_indexes(j,i) = min(find(diffTheta==min(diffTheta)));
									end
								end
							end
						end
						
						for i=tplot % build the S matrix for plotting
							for j=1:n
								ti=Thetalevel_indexes(i,j);
								if ~isnan(ti)
									interval=max(ti-1,1):min(ti+1,m); %interval is one above and one below ti
									a = trimPTMP(ti,j) - trimPTMP(interval, j);
									if( trimPTMP(ti,j)>tlevels(i) )
										gg=find(a>0);
										if( ~isempty(gg) )
											b=find(a==min(a(gg))); %find the level with min +ve diff
											ki=interval(b);
										else
											ki=ti;
										end
									end
									if( trimPTMP(ti,j)<tlevels(i) )
										gg=find(a<0);
										if( ~isempty(gg) )
											b=find(-a==min(-a(gg))); %find the level with min -ve diff
											ki=interval(b);
										else
											ki=ti;
										end
									end
									if( trimPTMP(ti,j)==tlevels(i) )
										ki=ti;
									end
									if( ki~=ti&~isnan(trimSAL(ti,j))&~isnan(trimSAL(ki,j))&~isnan(trimPTMP(ti,j))&~isnan(trimPTMP(ki,j)) )
										Sint(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trimSAL(ti,j), trimSAL(ki,j)], tlevels(i) );
									else
										Sint(i,j) = trimSAL(ti,j); % interpolate if possible because that is more accurate than using closest points
									end
									if( ki~=ti&~isnan(trim_mapped_sal(ti,j))&~isnan(trim_mapped_sal(ki,j))&~isnan(trimPTMP(ti,j))&~isnan(trimPTMP(ki,j)) )
										Smap(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_mapped_sal(ti,j), trim_mapped_sal(ki,j)], tlevels(i) );
										Smaperr(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_mapsalerrors(ti,j), trim_mapsalerrors(ki,j)], tlevels(i) );
									else
										Smap(i,j)=trim_mapped_sal(ti,j); % interpolate if possible because that is more accurate than using closest points
										Smaperr(i,j)=trim_mapsalerrors(ti,j); % interpolate if possible because that is more accurate than using closest points
									end
									if( ki~=ti&~isnan(trim_cal_SAL(ti,j))&~isnan(trim_cal_SAL(ki,j))&~isnan(trimPTMP(ti,j))&~isnan(trimPTMP(ki,j)) )
										Scal(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_cal_SAL(ti,j), trim_cal_SAL(ki,j)], tlevels(i) );
										Scalerr(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_cal_SAL_err(ti,j), trim_cal_SAL_err(ki,j)], tlevels(i) );
									else
										Scal(i,j)=trim_cal_SAL(ti,j); % interpolate if possible because that is more accurate than using closest points
										Scalerr(i,j)=trim_cal_SAL_err(ti,j); % interpolate if possible because that is more accurate than using closest points
									end
								end
							end
						end
					end
				end
						
				for j=1:n
					ii=find(~isnan(Smaperr(:,j)));
					cov = build_ptmp_cov( tlevels(ii));
					mean_err_mapp(j)=sum(sum(diag(Smaperr(ii,j))*cov))/(length(ii)*length(ii));
				end
				
				
				
				figure(2)
				 
				%keyboard
				%subplot(2,1,1)
				%plot(PROFILE_NO, avg_Soffset, 'b-');
				hold on
				%plot(PROFILE_NO, avg_Soffset, 'g-');
				% Plot station by station fit
				ok = find(isfinite(avg_Staoffset));
				 %keyboard
				addpath('/home1/homedir5/perso/ccabanes/matlab/lib_downlaod/')
				coly=[1 0 0];
				colerr=[1 0.8 0.8];
				%colerr=[0.8 0.8 1];
				 if THEPLOT==1
				 coly=[1 0 0];
				colerr=[1 0.8 0.8];
				 myebcc(PROFILE_NO(ok),avg_Staoffset(ok),mean_err_mapp(ok),coly,colerr)
				elseif THEPLOT==2
				coly=[0 0 1];
				colerr=[0.8 0.8 1];
				myebcc(PROFILE_NO(ok),avg_Staoffset(ok),mean_err_mapp(ok),coly,colerr)
				 % % plot(PROFILE_NO(ok), -avg_Staoffset(ok)+mean_err_mapp(ok), 'b-','LineWidth',0.5); 
				 % % plot(PROFILE_NO(ok), -avg_Staoffset(ok)-mean_err_mapp(ok), 'b-','LineWidth',0.5); 
				 end

				
				
				if THEPLOT==1
				profile1=PROFILE_NO(ok);
				offset1=avg_Staoffset(ok);
				%p1=plot(PROFILE_NO(ok), -avg_Staoffset(ok), 'r-','LineWidth',1.5);
				elseif THEPLOT==2
				profile2=PROFILE_NO(ok);
				offset2=avg_Staoffset(ok);
				%p2=plot(PROFILE_NO(ok), -avg_Staoffset(ok), 'b-','LineWidth',2);
				end
				
				%legend('2 x cal error','1 x cal error',['1-1 profile fit (' strrep(po_system_configuration.CONFIGURATION_FILE(1:end-4),'_',' ') ')' ], 'Old DM correction', 'Location', 'SouthOutside','orientation','horizontal');
				
			   % errorbar(PROFILE_NO, avg_Soffset, 2*avg_Soffset_err,'b')
			   % errorbar(PROFILE_NO, avg_Soffset, avg_Soffset_err,'go-')
				%plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r-');
				
				% %keyboard
				% addpath('/home1/homedir5/perso/ccabanes/matlab/lib_downlaod/')
				% coly=[1 0 0];
				% colerr=[1 0.8 0.8];
				% %colerr=[0.8 0.8 1];
				 % if THEPLOT==1
				 % coly=[1 0 0];
				% colerr=[1 0.8 0.8];
				 % myebcc(PROFILE_NO(ok),-avg_Staoffset(ok),mean_err_mapp(ok),coly,colerr)
				% elseif THEPLOT==2
				% coly=[0 0 1];
				% colerr=[0.8 0.8 1];
				% myebcc(PROFILE_NO(ok),-avg_Staoffset(ok),mean_err_mapp(ok),coly,colerr)
				 % % % plot(PROFILE_NO(ok), -avg_Staoffset(ok)+mean_err_mapp(ok), 'b-','LineWidth',0.5); 
				 % % % plot(PROFILE_NO(ok), -avg_Staoffset(ok)-mean_err_mapp(ok), 'b-','LineWidth',0.5); 
				 % end
				%axis([ 0, max(PROFILE_NO)+1, min([avg_Soffset-avg_Soffset_err,0])-.02, max([avg_Soffset+avg_Soffset_err,0])+.02 ])
				%axis([ 0, max(PROFILE_NO)+1, min(avg_Staoffset), max(avg_Staoffset)])
				
				%plot(PROFILE_NO, avg_Soffset, 'b-');
				%hold on
				%plot(PROFILE_NO, avg_Soffset, 'g-');
				%errorbar(PROFILE_NO, avg_Soffset, 2*avg_Soffset_err,'b')
				%errorbar(PROFILE_NO, avg_Soffset, avg_Soffset_err,'go-')
				
				axis tight
				plot( [0, max(PROFILE_NO)+1], [0,0], 'k-')
				set(gca,'FontSize',10)
				%xlabel(['Float profile number (date of the 1^{st} cycle: ' thedate_1cycle ')' ]);
				xlabel(['Float profile number ' ]);
				ylabel('\Delta S(PSS-78) ')
				%title( [' vertically-averaged salinity (PSS-78) additive correction \Delta S_{sta} (red: config', strrep(po_system_configuration.CONFIGURATION_FILE(1:end-4),'_',' '), ')'] );
				%title( ['B) Vertically-averaged salinity  additive correction and model fit'])
				title( ['cycle-by-cycle model drift'])
				grid on
				box on
				%plot(F.cycle_number.data,mean_diff,'m-','LineWidth',2)
				
				 % hold on
				 % if THEPLOT==1
				% %     plot(PROFILE_NO(ok), -avg_Staoffset(ok)-mean_err_mapp(ok), 'r-','LineWidth',0.5);
				% %     plot(PROFILE_NO(ok), -avg_Staoffset(ok)+mean_err_mapp(ok), 'r-','LineWidth',0.5);				 
				 % elseif THEPLOT==2
				 % plot(PROFILE_NO(ok), avg_Staoffset(ok)-mean_err_mapp(ok), 'b-','LineWidth',1); 
				 % plot(PROFILE_NO(ok), avg_Staoffset(ok)+mean_err_mapp(ok), 'b-','LineWidth',1); 
				 % end
				 grid on
				 box on
			end %if(isempty(find(isnan(cal_SAL)==0))==0) ---------------
			
			
			%drawnow
			%set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
			%print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_9.eps'));
			%%print('-dpdf ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_9.pdf'));
			
			%  set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
			%  print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_9.eps'));
			%  disp(['!ps2pdf ' po_system_configuration.FLOAT_PLOTS_DIRECTORY pn_float_dir pn_float_name '_9.eps'])
			%  eval(['!ps2pdf ' po_system_configuration.FLOAT_PLOTS_DIRECTORY pn_float_dir pn_float_name '_9.eps ' po_system_configuration.FLOAT_PLOTS_DIRECTORY pn_float_dir pn_float_name '_9.pdf'  ])
			
			
			
			
			
		end %if(isempty(find(isnan(mapped_sal)==0))==0) -----------
		
	end %if(isempty(find(isnan(PRES)==0))==0) ------------------
end

p1=plot(profile1,offset1, 'r-','LineWidth',2);
p2=plot(profile2,offset2, 'b-','LineWidth',2);

h=[p1,p2];
legend(h,{strrep(POSTFIX_D(2:end),'_','-'),strrep(POSTFIX_S(2:end),'_','-')})
%legend(h,{'3000-4000','400-800'})

file_out=[C.DIR_DATA 'float_plots/CONFIG' num2str(configow) '/' tabfloat '/' tabfloat '_15.pdf'];
print ('-dpdf' ,file_out)
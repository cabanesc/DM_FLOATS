
function plot_diagnostics_ow_figure9_4( pn_float_dir, pn_float_name, po_system_configuration,dacname,C )

%
% Annie Wong, 14 June 2011
% Breck Owens, October 2006
%--------------------------------------------------------------------------

%pn_float_dir='uw/';
%pn_float_name='R5902134';
%po_system_configuration = load_configuration( 'ow_config.txt' );


%close all

DIR_FTP=C.DIR_FTP;


% Move in ini_path.m by T. Reynaud 21.09.2020
%addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/Lib/ObsInSitu/General/')
%addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/Lib/ObsInSitu/General/Coriolis/')
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
            strcat( po_system_configuration.FLOAT_CALSERIES_PREFIX , pn_float_name, po_system_configuration.FLOAT_MAPPED_POSTFIX ) ) );
        
        use_theta_gt = lo_float_calseries.use_theta_gt;
        use_theta_lt = lo_float_calseries.use_theta_lt;
        use_pres_gt = lo_float_calseries.use_pres_gt;
        use_pres_lt = lo_float_calseries.use_pres_lt;
        use_percent_gt = lo_float_calseries.use_percent_gt;
        
        % load the station by station fits
        load(fullfile( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir,...
            strcat( po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ),'-regexp','^sta')
        
        
        % plot the float locations (figure 1) -----------------------
        
        load( fullfile( po_system_configuration.CONFIG_DIRECTORY, po_system_configuration.CONFIG_COASTLINES ), 'coastdata_x', 'coastdata_y' );
        
        
        %keyboard
        [m,n] = size(PRES);
        
        figure
        subplot(2,1,1)
        
        
        %Topo_ficin='/home/lpoargo1/DMARGO/OW/TOPO/topo.onetenthdeg.atl.nc';
        %Topo_ficin='/Users/thierry_reynaud/IFREMER/MATLAB/DM_FLOATS/CFLAG_EXTERNAL_DATA/TOPOGRAPHY/topo.onetenthdeg.nc';
        Topo_ficin=C.FILE_TOPO;
        Topo=read_netcdf_allthefile(Topo_ficin);
        hold on
        
        %load colormap.mat
        
        
        grid on
        %colormap(w)
        box on
        
        % %  set(gcf,'defaultaxeslinewidth',2)
        % %  set(gcf,'defaultlinelinewidth',2)
        % %  set(gcf,'defaultaxesfontsize',12)
        
        colormap(jet(n));
        c=colormap;
        
        x=[];
        y=[];
        if(isempty(selected_hist)==0)
            HIST.longitude.data=selected_hist(:,1);
            HIST.latitude.data=selected_hist(:,2);
            HIST=shiftEW(HIST,'longitude','grwch');
            x=HIST.longitude.data;
            y=HIST.latitude.data;
        end
        
        plot(x,y,'.','Color',[0.6 0.6 0.6],'MarkerSize',8)
        hold on
        plot(LONG,LAT,'r-','LineWidth',1);
        
        
        
        minx = max(min(min(x),min(LONG))-5,-100);
        maxx = max(max(x),max(LONG))+5;
        miny = min(min(y),min(LAT))-5;
        maxy=  max(max(y),max(LAT))+5;
        
        if maxy<75; maxy=maxy;end
        
        if  maxy>=75; maxy=75;end
        
        
        axis([minx maxx miny maxy])
        
        %legend('float','historical points','Location','Best')
        %plot(coastdata_x,coastdata_y,'k.-');
        
        for i=1:n
            h=plot(LONG(i),LAT(i),'+','MarkerSize',5);
            set(h,'color',c(i,:));
            % j=text(LONG(i),LAT(i),int2str(PROFILE_NO(i)));
            %set(j,'color',c(i,:),'fontsize',12,'hor','cen');
        end
        
        
        set(gca,'FontSize',11)
        xlabel('Longitude');
        ylabel('Latitude');
        
        
        
        %title( strcat(title_floatname, ' ARGO +CTD REF : ', po_system_configuration.CONFIGURATION_FILE(1:end-4), ',  (', strtrim(F.pi_name.data(1,:)), '---', strtrim(F.data_centre.data(1,:)),')' ),'interpreter','none' );
        
        %title( ['A) Argo float ' title_floatname ': path and reference data (ARGO + CTD)']);
        %keyboard
        %title( ['A) Argo float ' title_floatname ', ' strtrim(F.pi_name.data(1,:)) ', ' strtrim(F.data_centre.data(1,:))]);
        iwrap=1;
        
        lat=Topo.lat.data;
        bathy.topo.data=NaN*ones(size(Topo.topo.data));
        if iwrap==1
            idx_ga=find(Topo.lon.data>=180);
            idx_dr=find(Topo.lon.data<180);
            if ~isempty(idx_ga)
                lon=NaN*ones(size(Topo.lon.data));
                lon(1:length(idx_ga))=Topo.lon.data(idx_ga)-360;
                lon(length(idx_ga)+1:end)=Topo.lon.data(idx_dr);
                
                bathy.topo.data(:,1:length(idx_ga))=Topo.topo.data(:,idx_ga)-360;
                bathy.topo.data(:,length(idx_ga)+1:end)=Topo.topo.data(:,idx_dr);
                
            end
        else
            bathy.topo.data=Topo.topo.data;
            lon=Topo.lon.data;
        end
        [bathy.lon.data,bathy.lat.data]=meshgrid(lon,lat);
        
        iix=find(bathy.lon.data(1,:)>minx & bathy.lon.data(1,:)<maxx);
        jjy=find(bathy.lat.data(:,1)>miny & bathy.lat.data(:,1)<maxy);
        
        
        
        % Commented by T. Reynaud 02/11/2020
        %         iix=find(Topo.lon.data(1,:)>minx&Topo.lon.data(1,:)<maxx);
        %         jjy=find(Topo.lat.data(:,1)>miny&Topo.lat.data(:,1)<maxy);
        %         contour(Topo.lon.data(jjy,iix),Topo.lat.data(jjy,iix),Topo.topo.data(jjy,iix),[-1000 -1000],'LineColor',[0.6 0.5 0.4]);
        %         contour(Topo.lon.data(jjy,iix),Topo.lat.data(jjy,iix),Topo.topo.data(jjy,iix),[-2000 -2000],'LineColor',[0.6 0.5 0.4]);
        %         contour(Topo.lon.data(jjy,iix),Topo.lat.data(jjy,iix),Topo.topo.data(jjy,iix),[-3000 -3000],'LineColor',[0.8 0.7 0.6]);
        %         contour(Topo.lon.data(jjy,iix),Topo.lat.data(jjy,iix),Topo.topo.data(jjy,iix),[-4000 -4000],'LineColor',[0.8 0.8 0.8]);
        %         contour(Topo.lon.data(jjy,iix),Topo.lat.data(jjy,iix),Topo.topo.data(jjy,iix),[0 0],'LineColor',[0 0 0]);
        
        % Added by T. Reynaud 02/11/2020
        contour(bathy.lon.data(jjy,iix),bathy.lat.data(jjy,iix),bathy.topo.data(jjy,iix),[-1000 -1000],'LineColor',[0.6 0.5 0.4]);
        contour(bathy.lon.data(jjy,iix),bathy.lat.data(jjy,iix),bathy.topo.data(jjy,iix),[-2000 -2000],'LineColor',[0.6 0.5 0.4]);
        contour(bathy.lon.data(jjy,iix),bathy.lat.data(jjy,iix),bathy.topo.data(jjy,iix),[-3000 -3000],'LineColor',[0.8 0.7 0.6]);
        contour(bathy.lon.data(jjy,iix),bathy.lat.data(jjy,iix),bathy.topo.data(jjy,iix),[-4000 -4000],'LineColor',[0.8 0.8 0.8]);
        contour(bathy.lon.data(jjy,iix),bathy.lat.data(jjy,iix),bathy.topo.data(jjy,iix),[0 0],'LineColor',[0 0 0]);
        
        
        
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
            
            
            
            
            
            subplot(2,1,2)
            plot(PROFILE_NO, avg_Soffset, 'b-');
            hold on
            plot(PROFILE_NO, avg_Soffset, 'g-');
            % Plot station by station fit
            ok = find(isfinite(avg_Staoffset));
            plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r-+','LineWidth',1.5);
            
            
            %legend('2 x cal error','1 x cal error',['1-1 profile fit (' strrep(po_system_configuration.CONFIGURATION_FILE(1:end-4),'_',' ') ')' ], 'Old DM correction', 'Location', 'SouthOutside','orientation','horizontal');
            
            errorbar(PROFILE_NO, avg_Soffset, 2*avg_Soffset_err,'b')
            errorbar(PROFILE_NO, avg_Soffset, avg_Soffset_err,'go-')
            %plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r-');
            
            %keyboard
            %addpath('/home1/homedir5/perso/ccabanes/matlab/lib_downlaod/')
            coly=[1 0 0];
            colerr=[1 0.8 0.8];
            
            myebcc(PROFILE_NO(ok),avg_Staoffset(ok),mean_err_mapp(ok),coly,colerr)
            
            %axis([ 0, max(PROFILE_NO)+1, min([avg_Soffset-avg_Soffset_err,0])-.02, max([avg_Soffset+avg_Soffset_err,0])+.02 ])
            %axis([ 0, max(PROFILE_NO)+1, min(avg_Staoffset), max(avg_Staoffset)])
            
            plot(PROFILE_NO, avg_Soffset, 'b-');
            hold on
            plot(PROFILE_NO, avg_Soffset, 'g-');
            errorbar(PROFILE_NO, avg_Soffset, 2*avg_Soffset_err,'b')
            errorbar(PROFILE_NO, avg_Soffset, avg_Soffset_err,'go-')
            
            axis tight
            plot( [0, max(PROFILE_NO)+1], [0,0], 'k-')
            set(gca,'FontSize',10)
            %xlabel(['Float profile number (date of the 1^{st} cycle: ' thedate_1cycle ')' ]);
            xlabel(['Float profile number ' ]);
            ylabel('\Delta S(PSS-78) ')
            %title( [' vertically-averaged salinity (PSS-78) additive correction \Delta S_{sta} (red: config', strrep(po_system_configuration.CONFIGURATION_FILE(1:end-4),'_',' '), ')'] );
            title( ['B) Vertically-averaged salinity  additive correction and model fit'])
            grid on
            
            
            if C.plotDm==1
                filename = [DIR_FTP strtrim(dacname) '/' strtrim(title_floatname) '/' strtrim(title_floatname) '_prof.nc'];
                % disp([num2str(unum(k)) ' ' F.wmo_inst_type.data(1,:)])
                
                F = read_netcdf_allthefile(filename);
                F = replace_fill_bynan(F);
                F = format_flags_char2num(F);
                param='psal';
                paramad='psal_adjusted';
                paramqc='psal_qc';
                paramadqc='psal_adjusted_qc';

                N_PROF_this_file = size(F.cycle_number.data,1);

                mean_diff = NaN*zeros(N_PROF_this_file,1);
                
                %isnanprof=isnan(F.(param).data)|isnan(F.(paramad).data)|F.(param).data>9999999|F.(paramad).data>9999999;
                isnanprof=isnan(F.(param).data)|isnan(F.(paramad).data);
                F.(paramad).data(isnanprof)=NaN;
                F.(param).data(isnanprof)=NaN;
                
                %F.(param).data(F.(param).data>9999999)=NaN;
                %F.(paramad).data(F.(paramad).data>9999999)=NaN;
                
                mean_diff = meanoutnan(F.(paramad).data,2)-meanoutnan(F.(param).data,2);
                mean_diff = meanoutnan(F.(paramad).data,2)-meanoutnan(F.(param).data,2);
                plot_cycle=[min(F.cycle_number.data):max(F.cycle_number.data)];
                plot_mean=NaN*ones(size(plot_cycle));
                [iiu,iiy]=ismember(F.cycle_number.data,plot_cycle);
                plot_mean(iiy)=mean_diff;
                plot(plot_cycle,plot_mean,'m-','LineWidth',2)
            end
            
        end %if(isempty(find(isnan(cal_SAL)==0))==0) ---------------
        
        
        drawnow
        set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
        print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_9.eps'));
        print('-dpdf', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_9.pdf'));
        
        %  set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
        %  print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_9.eps'));
        %  disp(['!ps2pdf ' po_system_configuration.FLOAT_PLOTS_DIRECTORY pn_float_dir pn_float_name '_9.eps'])
        %  eval(['!ps2pdf ' po_system_configuration.FLOAT_PLOTS_DIRECTORY pn_float_dir pn_float_name '_9.eps ' po_system_configuration.FLOAT_PLOTS_DIRECTORY pn_float_dir pn_float_name '_9.pdf'  ])
        
        
        
        disp(' ')
        disp(' ')
        tlevels
        filename_stat=fullfile( po_system_configuration.FLOAT_MAPPED_DIRECTORY, pn_float_dir, ...
            strcat( pn_float_name, '_stat', po_system_configuration.FLOAT_MAPPED_POSTFIX ) );
        if exist(filename_stat)
            lo_float_source_stat_data = load(filename_stat);
            ii_large_all=NaN*zeros(length(PROFILE_NO),1);
            ii_small_all=NaN*zeros(length(PROFILE_NO),1);
            ii_large_5=NaN*zeros(length(PROFILE_NO),1);
            ii_large_10=NaN*zeros(length(PROFILE_NO),1);
            ii_large_2=NaN*zeros(length(PROFILE_NO),1);
            ii_large_1=NaN*zeros(length(PROFILE_NO),1);
            ii_small_10=NaN*zeros(length(PROFILE_NO),1);
            ii_small_5=NaN*zeros(length(PROFILE_NO),1);
            ii_small_2=NaN*zeros(length(PROFILE_NO),1);
            ii_small_1=NaN*zeros(length(PROFILE_NO),1);
            for i=1:length(PROFILE_NO)
                
                iprof=lo_float_source_stat_data.profile_index(i);
                if ~isnan(iprof)
                    itruc=ismember(lo_float_source_stat_data.index_un,lo_float_source_stat_data.index_each{i});
                    
                    la_bhist_lat = lo_float_source_stat_data.la_bhist_lat_all(itruc);
                    la_bhist_long = lo_float_source_stat_data.la_bhist_long_all(itruc);
                    la_bhist_dates = lo_float_source_stat_data.la_bhist_dates_all(itruc);
                    la_bhist_Z = lo_float_source_stat_data.la_bhist_Z_all(itruc);
                    Min_ptmp = lo_float_source_stat_data.Min_ptmp(itruc)';
                    
                    fl_LAT = lo_float_source_data.LAT( iprof) ;
                    fl_LONG = lo_float_source_data.LONG( iprof) ;
                    fl_DATES = lo_float_source_data.DATES(iprof) ;
                    
                    if(fl_LONG>180) % m_tbase inputs longitudes from 0 to +/- 180
                        fl_LONG1=fl_LONG-360;
                    else
                        fl_LONG1=fl_LONG;
                    end
                    if(la_bhist_long>180) % m_tbase inputs longitudes from 0 to +/- 180
                        la_bhist_long1=la_bhist_long-360;
                    else
                        la_bhist_long1=la_bhist_long;
                    end
                    
                    m_proj('mercator','long', [min(fl_LONG1)-1, max(fl_LONG1)+1], 'lat', [min(fl_LAT)-1, max(fl_LAT)+1] );
                    [elev,x,y] = m_tbase( [min(fl_LONG1)-1, max(fl_LONG1)+1, min(fl_LAT)-1, max(fl_LAT)+1] );
                    fl_Z = -interp2( x,y,elev, fl_LONG1, fl_LAT, 'linear'); % -ve bathy values
                    
                    PV_float = (2*7.292*10^-5.*sin(fl_LAT.*pi/180))./fl_Z;
                    PV_bhist = (2*7.292*10^-5.*sin(la_bhist_lat.*pi/180))./la_bhist_Z;
                    
                    if( str2num(po_system_configuration.MAP_USE_PV)==1 ) % if PV is wanted
                        ellipse_large = sqrt( (la_bhist_long1-fl_LONG1).^2./(str2num(po_system_configuration.MAPSCALE_LONGITUDE_LARGE)).^2 + (la_bhist_lat-fl_LAT).^2./(str2num(po_system_configuration.MAPSCALE_LATITUDE_LARGE)).^2 +...
                            ((PV_float-PV_bhist)./sqrt( PV_float.^2+PV_bhist.^2 )./str2num(po_system_configuration.MAPSCALE_PHI_LARGE)).^2 ) ;
                        ellipse_small = sqrt( (la_bhist_long1-fl_LONG1).^2./(str2num(po_system_configuration.MAPSCALE_LONGITUDE_SMALL)).^2 + (la_bhist_lat-fl_LAT).^2./(str2num(po_system_configuration.MAPSCALE_LATITUDE_SMALL)).^2 +...
                            ((PV_float-PV_bhist)./sqrt( PV_float.^2+PV_bhist.^2 )./str2num(po_system_configuration.MAPSCALE_PHI_SMALL)).^2 ) ;
                    else % if PV is unwanted ---
                        ellipse_large = sqrt( (la_bhist_long1-fl_LONG1).^2./(str2num(po_system_configuration.MAPSCALE_LONGITUDE_LARGE)).^2 + (la_bhist_lat-fl_LAT).^2./(str2num(po_system_configuration.MAPSCALE_LATITUDE_LARGE)).^2 );
                        ellipse_small = sqrt( (la_bhist_long1-fl_LONG1).^2./(str2num(po_system_configuration.MAPSCALE_LONGITUDE_SMALL)).^2 + (la_bhist_lat-fl_LAT).^2./(str2num(po_system_configuration.MAPSCALE_LATITUDE_SMALL)).^2 );
                    end
                    
                    ii_large_1(i)=sum(ellipse_large<1&abs(la_bhist_dates-fl_DATES)<1&Min_ptmp<medianoutnan(tlevels));
                    ii_large_2(i)=sum(ellipse_large<1&abs(la_bhist_dates-fl_DATES)<2&Min_ptmp<medianoutnan(tlevels));
                    ii_large_5(i)=sum(ellipse_large<1&abs(la_bhist_dates-fl_DATES)<5&Min_ptmp<medianoutnan(tlevels));
                    ii_large_10(i)=sum(ellipse_large<1&abs(la_bhist_dates-fl_DATES)<10&Min_ptmp<medianoutnan(tlevels));
                    ii_large_all(i)=sum(ellipse_large<1&Min_ptmp<medianoutnan(tlevels));
                    
                    ii_small_all(i)=sum(ellipse_small<1&Min_ptmp<medianoutnan(tlevels));
                    ii_small_10(i)=sum(ellipse_small<1 & abs(la_bhist_dates-fl_DATES)<10&Min_ptmp<medianoutnan(tlevels));
                    ii_small_5(i)=sum(ellipse_small<1 & abs(la_bhist_dates-fl_DATES)<5&Min_ptmp<medianoutnan(tlevels));
                    ii_small_2(i)=sum(ellipse_small<1 & abs(la_bhist_dates-fl_DATES)<2&Min_ptmp<medianoutnan(tlevels));
                    ii_small_1(i)=sum(ellipse_small<1 & abs(la_bhist_dates-fl_DATES)<1&Min_ptmp<medianoutnan(tlevels));
                end
            end
        end
        ii_large=[ii_large_1,ii_large_2,ii_large_5,ii_large_10,ii_large_all,ii_large_all];
        ii_small=[ii_small_1,ii_small_2,ii_small_5,ii_small_10,ii_small_all,ii_small_all];
        ii_large_nan=ii_large;
        ii_large_nan(ii_large==0)=NaN;
        ii_small_nan=ii_small;
        ii_small_nan(ii_small==0)=NaN;
        
        fig=figure;
        
        
        %subplot (2,15,[3:15])
        subplot ('Position',[0.22 0.52 0.7 0.35])
        %keyboard
        pcolor(PROFILE_NO,[1:6],ii_large_nan');
        %sanePColor(PROFILE_NO,[1:5],ii_large_nan');
        set(gca,'YtickLabel',{''})
        set(gca,'XtickLabel',{''})
        box on
        colormap(c)
        title('within LARGE spatial scales ')
        %xlabel('Profile Number')
        colorbar
        hold on
        box on
        %subplot (2,15,[1:2])
        set(gca,'TickDir','out')
        
        %subplot ('Position',[0.07 0.52 0.1 0.35])
        subplot ('Position',[0.1 0.52 0.1 0.35])
        box on
        hold on
        grid on
        pcolor([meanoutnan(ii_large);meanoutnan(ii_large)]')
        %sanePColor([meanoutnan(ii_large);meanoutnan(ii_large)]')
        caxis([0 max(max(ii_large))])
        c=cool;
        
        %axis([1 2 0.5  5.5])
        %axis([0.5 2.5 0.5  5.5])
        set(gca,'Ytick',[1.5 2.5 3.5 4.5 5.5])
        set(gca,'YtickLabel',[{'1','2','5','10','all',''}])
        
        ylabel('TIME SCALE ( yr )')
        %xlabel('average')
        for i=1:5
            text(1.5,i+0.5,[num2str(floor(meanoutnan(ii_large(:,i)))) ' +/- ' num2str(floor(stdoutnan(ii_large(:,i))))] ,'color','w','fontweight','bold','fontsize',7,'HorizontalAlignment','center')
        end
        
        %subplot (2,15,[18:30])
        subplot ('Position',[0.22 0.1 0.7 0.35])
        pcolor(PROFILE_NO,[1:6],ii_small_nan');
        %sanePColor(PROFILE_NO,[1:5],ii_small_nan');
        set(gca,'YtickLabel',{''})
        box on
        colormap(c)
        title(' within SMALL spatial scales ')
        xlabel('Profile Number')
        colorbar
        hold on
        box on
        set(gca,'TickDir','out')
        
        %subplot (2,15,[16:17])
        subplot ('Position',[0.1 0.1 0.1 0.35])
        box on
        hold on
        grid on
        pcolor([meanoutnan(ii_small);meanoutnan(ii_small)]')
        %sanePColor([meanoutnan(ii_small);meanoutnan(ii_small)]')
        caxis([0 max(max(ii_small))])
        %axis([1 2 0.5  5.5])
        %axis([0.5 2.5 0.5  5.5])
        set(gca,'Ytick',[1.5 2.5 3.5 4.5 5.5])
        set(gca,'YtickLabel',[{'1','2','5','10','all',''}])
        ylabel('TIME SCALE ( yr )')
        xlabel('Average (+/- std)')
        for i=1:5
            text(1.5,i+0.5,[num2str(floor(meanoutnan(ii_small(:,i)))) ' +/- ' num2str(floor(stdoutnan(ii_small(:,i))))] ,'color','w','fontweight','bold','fontsize',7,'HorizontalAlignment','center')
        end
        
        ax=axes('Units','Normal','Position',[.01 .01 .95 0.90],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        title(['Number of reference data available for float ' pn_float_name ' at theta = ' num2str(floor(medianoutnan(tlevels)*100)/100)],'FontWeight','bold','FontSize',10)
        
        drawnow
        set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,7,5]);
        print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_10.eps'));
        print('-dpng', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_10.png'));
        print('-dpdf', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_10.pdf'));
        
        disp(' ')
        disp(' ')
        tlevels
        
        filename_stat=fullfile( po_system_configuration.FLOAT_MAPPED_DIRECTORY, pn_float_dir, ...
            strcat( pn_float_name, '_stat', po_system_configuration.FLOAT_MAPPED_POSTFIX ) );
        if exist(filename_stat)
            lo_float_source_stat_data = load(filename_stat);
            ii_large_all=NaN*zeros(length(PROFILE_NO),1);
            ii_small_all=NaN*zeros(length(PROFILE_NO),1);
            ii_large_5=NaN*zeros(length(PROFILE_NO),1);
            ii_large_10=NaN*zeros(length(PROFILE_NO),1);
            ii_large_2=NaN*zeros(length(PROFILE_NO),1);
            ii_large_1=NaN*zeros(length(PROFILE_NO),1);
            ii_small_10=NaN*zeros(length(PROFILE_NO),1);
            ii_small_5=NaN*zeros(length(PROFILE_NO),1);
            ii_small_2=NaN*zeros(length(PROFILE_NO),1);
            ii_small_1=NaN*zeros(length(PROFILE_NO),1);
            for i=1:length(PROFILE_NO)
                
                iprof=lo_float_source_stat_data.profile_index(i);
                if ~isnan(iprof)
                    itruc=ismember(lo_float_source_stat_data.index_un,lo_float_source_stat_data.index_each{i});
                    
                    la_bhist_lat = lo_float_source_stat_data.la_bhist_lat_all(itruc);
                    la_bhist_long = lo_float_source_stat_data.la_bhist_long_all(itruc);
                    la_bhist_dates = lo_float_source_stat_data.la_bhist_dates_all(itruc);
                    la_bhist_Z = lo_float_source_stat_data.la_bhist_Z_all(itruc);
                    Min_ptmp = lo_float_source_stat_data.Min_ptmp(itruc)';
                    
                    fl_LAT = lo_float_source_data.LAT( iprof) ;
                    fl_LONG = lo_float_source_data.LONG( iprof) ;
                    fl_DATES = lo_float_source_data.DATES(iprof) ;
                    if(fl_LONG>180) % m_tbase inputs longitudes from 0 to +/- 180
                        fl_LONG1=fl_LONG-360;
                    else
                        fl_LONG1=fl_LONG;
                    end
                    if(la_bhist_long>180) % m_tbase inputs longitudes from 0 to +/- 180
                        la_bhist_long1=la_bhist_long-360;
                    else
                        la_bhist_long1=la_bhist_long;
                    end
                    m_proj('mercator','long', [min(fl_LONG1)-1, max(fl_LONG1)+1], 'lat', [min(fl_LAT)-1, max(fl_LAT)+1] );
                    [elev,x,y] = m_tbase( [min(fl_LONG1)-1, max(fl_LONG1)+1, min(fl_LAT)-1, max(fl_LAT)+1] );
                    fl_Z = -interp2( x,y,elev, fl_LONG1, fl_LAT, 'linear'); % -ve bathy values
                    
                    PV_float = (2*7.292*10^-5.*sin(fl_LAT.*pi/180))./fl_Z;
                    PV_bhist = (2*7.292*10^-5.*sin(la_bhist_lat.*pi/180))./la_bhist_Z;
                    
                    if( str2num(po_system_configuration.MAP_USE_PV)==1 ) % if PV is wanted
                        ellipse_large = sqrt( (la_bhist_long1-fl_LONG1).^2./(str2num(po_system_configuration.MAPSCALE_LONGITUDE_LARGE)).^2 + (la_bhist_lat-fl_LAT).^2./(str2num(po_system_configuration.MAPSCALE_LATITUDE_LARGE)).^2 +...
                            ((PV_float-PV_bhist)./sqrt( PV_float.^2+PV_bhist.^2 )./str2num(po_system_configuration.MAPSCALE_PHI_LARGE)).^2 ) ;
                        ellipse_small = sqrt( (la_bhist_long1-fl_LONG1).^2./(str2num(po_system_configuration.MAPSCALE_LONGITUDE_SMALL)).^2 + (la_bhist_lat-fl_LAT).^2./(str2num(po_system_configuration.MAPSCALE_LATITUDE_SMALL)).^2 +...
                            ((PV_float-PV_bhist)./sqrt( PV_float.^2+PV_bhist.^2 )./str2num(po_system_configuration.MAPSCALE_PHI_SMALL)).^2 ) ;
                    else % if PV is unwanted ---
                       ellipse_large = sqrt( (la_bhist_long1-fl_LONG1).^2./(str2num(po_system_configuration.MAPSCALE_LONGITUDE_LARGE)).^2 + (la_bhist_lat-fl_LAT).^2./(str2num(po_system_configuration.MAPSCALE_LATITUDE_LARGE)).^2 );
                        ellipse_small = sqrt( (la_bhist_long1-fl_LONG1).^2./(str2num(po_system_configuration.MAPSCALE_LONGITUDE_SMALL)).^2 + (la_bhist_lat-fl_LAT).^2./(str2num(po_system_configuration.MAPSCALE_LATITUDE_SMALL)).^2 );
                    
                        %ellipse = sqrt((grid_long-LONG).^2./(longitude_large*3).^2 + (grid_lat-LAT).^2./(latitude_large*3).^2) ;
                    end
                    
                    ii_large_1(i)=sum(ellipse_large<1&abs(la_bhist_dates-fl_DATES)<1&Min_ptmp<max(tlevels));
                    ii_large_2(i)=sum(ellipse_large<1&abs(la_bhist_dates-fl_DATES)<2&Min_ptmp<max(tlevels));
                    ii_large_5(i)=sum(ellipse_large<1&abs(la_bhist_dates-fl_DATES)<5&Min_ptmp<max(tlevels));
                    ii_large_10(i)=sum(ellipse_large<1&abs(la_bhist_dates-fl_DATES)<10&Min_ptmp<max(tlevels));
                    ii_large_all(i)=sum(ellipse_large<1&Min_ptmp<max(tlevels));
                    
                    ii_small_all(i)=sum(ellipse_small<1&Min_ptmp<max(tlevels));
                    ii_small_10(i)=sum(ellipse_small<1 & abs(la_bhist_dates-fl_DATES)<10&Min_ptmp<max(tlevels));
                    ii_small_5(i)=sum(ellipse_small<1 & abs(la_bhist_dates-fl_DATES)<5&Min_ptmp<max(tlevels));
                    ii_small_2(i)=sum(ellipse_small<1 & abs(la_bhist_dates-fl_DATES)<2&Min_ptmp<max(tlevels));
                    ii_small_1(i)=sum(ellipse_small<1 & abs(la_bhist_dates-fl_DATES)<1&Min_ptmp<max(tlevels));
                end
            end
        end
        ii_large=[ii_large_1,ii_large_2,ii_large_5,ii_large_10,ii_large_all,ii_large_all];
        ii_small=[ii_small_1,ii_small_2,ii_small_5,ii_small_10,ii_small_all,ii_small_all];
        ii_large_nan=ii_large;
        ii_large_nan(ii_large==0)=NaN;
        ii_small_nan=ii_small;
        ii_small_nan(ii_small==0)=NaN;
        
        fig=figure;
        
        
        %subplot (2,15,[3:15])
        subplot ('Position',[0.22 0.52 0.7 0.35])
        pcolor(PROFILE_NO,[1:6],ii_large_nan');
        %sanePColor(PROFILE_NO,[1:5],ii_large_nan');
        set(gca,'YtickLabel',{''})
        set(gca,'XtickLabel',{''})
        box on
        colormap(c)
        title('within LARGE spatial scales ')
        %xlabel('Profile Number')
        colorbar
        hold on
        box on
        %subplot (2,15,[1:2])
        set(gca,'TickDir','out')
        
        subplot ('Position',[0.1 0.52 0.1 0.35])
        box on
        hold on
        grid on
        pcolor([meanoutnan(ii_large);meanoutnan(ii_large)]')
        %sanePColor([meanoutnan(ii_large);meanoutnan(ii_large)]')
        caxis([0 max(max(ii_large))])
        c=cool;
        
        %axis([1 2 0.5  5.5])
        %axis([0.5 2.5 0.5  5.5])
        set(gca,'Ytick',[1.5 2.5 3.5 4.5 5.5])
        set(gca,'YtickLabel',[{'1','2','5','10','all',''}])
        
        ylabel('TIME SCALE ( yr )')
        %xlabel('average')
        for i=1:5
            text(1.5,i+0.5,[num2str(floor(meanoutnan(ii_large(:,i)))) ' +/- ' num2str(floor(stdoutnan(ii_large(:,i))))] ,'color','w','fontweight','bold','fontsize',7,'HorizontalAlignment','center')
        end
        
        %subplot (2,15,[18:30])
        subplot ('Position',[0.22 0.1 0.7 0.35])
        pcolor(PROFILE_NO,[1:6],ii_small_nan');
        %sanePColor(PROFILE_NO,[1:5],ii_small_nan');
        set(gca,'YtickLabel',{''})
        box on
        colormap(c)
        title(' within SMALL spatial scales ')
        xlabel('Profile Number')
        colorbar
        hold on
        box on
        set(gca,'TickDir','out')
        
        %subplot (2,15,[16:17])
        subplot ('Position',[0.1 0.1 0.1 0.35])
        box on
        hold on
        grid on
        pcolor([meanoutnan(ii_small);meanoutnan(ii_small)]')
        %sanePColor([meanoutnan(ii_small);meanoutnan(ii_small)]')
        caxis([0 max(max(ii_small))])
        %axis([1 2 0.5  5.5])
        %axis([0.5 2.5 0.5  5.5])
        set(gca,'Ytick',[1.5 2.5 3.5 4.5 5.5])
        set(gca,'YtickLabel',[{'1','2','5','10','all',''}])
        ylabel('TIME SCALE ( yr )')
        xlabel('Average (+/- std)')
        for i=1:5
            text(1.5,i+0.5,[num2str(floor(meanoutnan(ii_small(:,i)))) ' +/- ' num2str(floor(stdoutnan(ii_small(:,i))))] ,'color','w','fontweight','bold','fontsize',7,'HorizontalAlignment','center')
        end
        
        ax=axes('Units','Normal','Position',[.01 .01 .95 0.90],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        title(['Number of reference data available for float ' pn_float_name ' at theta = ' num2str(floor(max(tlevels)*100)/100)],'FontWeight','bold','FontSize',10)
        
        drawnow
        set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,7,5]);
        print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_11.eps'));
        print('-dpng', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_11.png'));
        print('-dpdf', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_11.pdf'));
    end
end %if(isempty(find(isnan(mapped_sal)==0))==0) -----------

end %if(isempty(find(isnan(PRES)==0))==0) ------------------



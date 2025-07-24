function [thetitle]=plot_traj_glob(FL,Topo,app)
% -========================================================
%   USAGE : plot_traj(FL)
%   PURPOSE : plot the trajectory of the float, on a map with bathy
% -----------------------------------
%   INPUT :
%     FL   (structure)  Float data
%   optional input :
%     Vec 1*n_cycle parameter to plot along the trajectory
%   OUTPUT :
%    thetitle (char) title of the plot
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx
%   CALLED SUBROUTINES: none
% ========================================================

if isfield(Topo,'lat')==0
    if isfield(Topo,'grd_lat')
        Topo.lat=Topo.grd_lat;
    end
end
if isfield(Topo,'lon')==0
    if isfield(Topo,'grd_lon')
        Topo.lon=Topo.grd_lon;
    end
end

Vec=[FL.min_cy:FL.max_cy];

% load bottom topography colormap
load ('./util/bathy_17_colormap.mat');

Topo = shiftEW(Topo,'lon','grwch');


if isfield(FL,'longitude')
    if isfield(FL,'latitude')
        
        
        FL = shiftEW(FL,'longitude','grwch');
        lonmin = floor(min(FL.longitude.data));
        lonmax = ceil(max(FL.longitude.data));
        latmin = floor(min(FL.latitude.data));
        latmax = ceil(max(FL.latitude.data));
        
        region.lonmin = max([lonmin-10,-180]);
        region.lonmax = min([lonmax+10,180]);
        region.latmin = max([latmin-10,-90]);
        region.latmax = min([latmax+10,90]);
        if region.latmin==70
            region.latmin=71
        end
        if lonmin<-100 && lonmax>100
            Topo = shiftEW(Topo,'lon','pacif');
            FL = shiftEW(FL,'longitude','pacif');
            lonmin = floor(min(FL.longitude.data));
            lonmax = ceil(max(FL.longitude.data));
            region.lonmin = max([lonmin-10,0]);
            region.lonmax = min([lonmax+10,360]);
        end
        
        %keyboard
        if min(size(Topo.lat.data))==1
            Topo.lat.data =repmat(Topo.lat.data,[1,length(Topo.lon.data)]);
            Topo.lon.data =repmat(Topo.lon.data',[size(Topo.lat.data,1),1]);
        end
        
        iy = (Topo.lat.data(:,1)>= region.latmin & Topo.lat.data(:,1)<=region.latmax );
        ix = (Topo.lon.data(1,:)>= region.lonmin & Topo.lon.data(1,:)<= region.lonmax);
        r=reshape(Topo.lat.data(iy,ix),sum(ix)*sum(iy),1);
        t=reshape(Topo.lon.data(iy,ix),sum(ix)*sum(iy),1);
        y=reshape(Topo.topo.data(iy,ix),sum(ix)*sum(iy),1);
        y(y>0)=NaN;
        
        if  latmin>=70|latmax>=80
            MAP_VISU.proj='stereographic';
            PARAM.isarctic=1;
            m_proj(MAP_VISU.proj,'lat',min((latmax+latmin)/2,90),'lon',max((lonmax+lonmin)/2,-180),'radius',min(latmax-latmin+5,30));
            [bathy,lon_bathy,lat_bathy]=m_tbase( [-180 180 max(latmin-10, -90) min(latmax+10,90)]);
            cvec=[-100000,-6500:500:1000,100000];
            [c,h]=m_contourf(lon_bathy,lat_bathy,bathy,cvec);
            %m_coast
            m_grid
        else
            MAP_VISU.proj ='miller'; MAP_VISU.reso='LR';
            PARAM.isarctic=0;
            m_proj(MAP_VISU.proj,'long',[region.lonmin,region.lonmax],'lati',[region.latmin,region.latmax]);
            cvec=[-100000,-6500:500:1000,100000];
            [c,h]=m_contourf(Topo.lon.data(iy,ix),Topo.lat.data(iy,ix),Topo.topo.data(iy,ix),cvec);
        end
        
        
        hold on;
                
        box on
        grid on
   
        
       
        set(h,'LineStyle','None');
        colormap(newmap);
%         p = get(h,'Children');
%         thechild=get(p,'CData');
%         cdat=cell2mat(thechild);
%         for i=1:length(cvec)-1
%             set(p(cdat>=cvec(i)& cdat< cvec(i+1)),'Facecolor',newmap(i,:),'LineStyle','none')
%         end
        %m_contour(Topo.lon.data(iy,ix),Topo.lat.data(iy,ix),Topo.topo.data(iy,ix),[0 0],'LineColor',[0.6 0.5 0.4]);
        
        %end
        
        Vec=Vec(~isnan(Vec));
        minvec=min(Vec);
        maxvec=max(Vec);
        col_pos=jet(maxvec-minvec+1);
        
        % Plot float trajectory
        m_plot(FL.longitude.data,FL.latitude.data,'b');
       
        %[poub,idx_cy,ipoub] = intersect(Vec,FL.cycle_number.data); % cc correction 19/10/2020
        [ispoub,idx_cy]=ismember(FL.cycle_number.data,Vec);
        if length(FL.cycle_number.data)<max(Vec)/3
            
            m_scatter(FL.longitude.data,FL.latitude.data,50,col_pos(idx_cy,:),'ok');
        end
        
        m_scatter(FL.longitude.data,FL.latitude.data,30,col_pos(idx_cy,:),'filled');
        
        if length(idx_cy)<max(Vec)/3
            m_plot(FL.longitude.data,FL.latitude.data,'w')
        end
        
        colorbar('delete')
        
        %c=colorbar('horiz')
        %xlabel(c,'Cycle Number')
        
        xth=get(gca,'XTick');
        yth=get(gca,'YTick');
        if PARAM.isarctic==0
            m_grid('xtick',4,'ytick',4);
            xlabel('Longitude')
            ylabel({'Latitude',' '})
        end
        thetitle = [' '];
        info1  = {'platform_number','pi_name','data_centre','juld'};
        info2  = {'Platform: ', ', PI: ', ', ',', date of the first cycle: '};
        
        %regionO=region;
        
        for k=1:length(info1)
            if isfield(FL,info1{k})
                if isfield(FL.(info1{k}), 'data')
                    if isequal(info1{k},'juld')
                        % transforme la date juld > jour dd/mm/yyyy
                        if ~isnan(FL.juld.data(1,:))
                        thedate = datestr( (FL.juld.data(1,:) + datenum('19500101','yyyymmdd')),'dd/mm/yyyy');
                        else
                         thedate ='NaN';
                        end
                        thetitle = [thetitle info2{k}  thedate  ];
                    else
                        thetitle = [thetitle info2{k}  strtrim(num2str(FL.(info1{k}).data(1,:)))  ];
                    end
                    
                end
            end
        end
        
        % title(thetitle,'interpreter','none')
        
    else
        error(' No latitude data in the float data structure')
    end
    
else
    error('No longitude data in the float data structure ')
end



function [thetitle]=plot_traj_glob_sans_bathy(FL,Vec)
% -========================================================
%   USAGE : plot_traj(FL)
%   PURPOSE : plot the trajectory of the float, on a map with bathy
% -----------------------------------
%   INPUT :
%     FL   (structure)  Float data
%   optional input :
%     Vec 1*n_cycle parameter to plot along the trajectory or color input in ('r','b','g','k','y','c','m')
%   OUTPUT :
%    thetitle (char) title of the plot
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx
%   CALLED SUBROUTINES: none
% ========================================================

%keyboard
 


if isfield(FL,'longitude')
    if isfield(FL,'latitude')
        
       
        %FL = shiftEW(FL,'longitude','grwch');
        
        region.lonmin = floor((min(FL.longitude.data)-10));
        region.lonmax = ceil(max(FL.longitude.data)+10);
        region.latmin = floor(min(FL.latitude.data)-10);
        region.latmax = ceil(max(FL.latitude.data)+10);
        
        
        
        
        
        
       
        
        %plot(FL.longitude.data,FL.latitude.data,'b')
        
        if  nargin==1
            Vec = FL.cycle_number.data;
        end
        %keyboard
        thedate = datevec(FL.juld.data+datenum('19500101','YYYYMMDD'));
        isnew=thedate(:,1)>=2010;
        %scatter(FL.longitude.data(isnew),FL.latitude.data(isnew),10,Vec(isnew),'filled')
        if ischar(Vec)
        scatter(FL.longitude.data,FL.latitude.data,30,Vec,'filled')
        elseif isnumeric(Vec)
        scatter(FL.longitude.data,FL.latitude.data,30,Vec,'filled')
        
        Vec=Vec(~isnan(Vec));
        minvec=min(Vec);
        maxvec=max(Vec);
        
        if minvec<maxvec
        caxis([minvec maxvec])
        else
        caxis([minvec maxvec+1])
        end
        end
        %xlabel('longitude')
        ylabel('latitude')
        theax=get(gca);
        theax.XAxisLocation='top';
        
        %colorbar
        
        
        thetitle = [' '];
        info1  = {'platform_number','pi_name','data_centre','juld'};
        info2  = {'Platform: ', ', PI: ', ', ',', date of the first cycle: '};
        
        %regionO=region;
        
        for k=1:length(info1)
            if isfield(FL,info1{k})
                if isfield(FL.(info1{k}), 'data')
                    if isequal(info1{k},'juld')
                        % transforme la date juld > jour dd/mm/yyyy
                        thedate = datestr( (FL.juld.data(1,:) + datenum('19500101','yyyymmdd')),'dd/mm/yyyy');
                        
                        thetitle = [thetitle info2{k}  thedate  ];
                    else
                        thetitle = [thetitle info2{k}  strtrim(num2str(FL.(info1{k}).data(1,:)))  ];
                    end
                    
                end
            end
        end
        
%        title(thetitle,'interpreter','none')
        
    else
        error(' No latitude data in the float data structure')
    end
    
else
    error('No longitude data in the float data structure ')
end




function [ pa_grid_lat, pa_grid_long, pa_grid_dates ] = get_region_ow(pa_wmo_numbers, pa_float_name, po_config_data) ;

% function [ pa_grid_lat, pa_grid_long, pa_grid_dates ] = get_region_ow( pa_wmo_numbers, pa_float_name, po_config_data ) ;
%
% This function gets historical data from the 25 selected WMO boxes,
% merges the CTD, BOT and Argo files, makes sure longitude is continuous
% around the 0-360 degree mark, and converts the dates of the
% historical data from time format2 to year.mo+ (changedates.m).
%
% pa_wmo_numbers can be NaN: when float profiles are out of range (65N, 65S),
% or when there's no .mat file in that box, denoted by 0 (e.g. box on land).
%
% Historical data have lat,long and dates organised in single rows.
% The output from this function gives lat,long and dates in columns,
% just for ease of checking .... really doesn't matter.
%
% Annie Wong, December 2007
% Breck Owens, December 2006
%keyboard

pa_grid_lat   = [ ] ;
pa_grid_long  = [ ] ;
pa_grid_dates = [ ] ;

[m,n]=size(pa_wmo_numbers);

% to save time: use of files that only contains lon and lat data
po_config_data.HISTORICAL_CTD_PREFIX='/historical_ctd/lonlat/ctd_';
po_config_data.HISTORICAL_ARGO_PREFIX='/argo_profiles/lonlat/argo_';
if exist([po_config_data.HISTORICAL_DIRECTORY '/historical_ctd/bad_ctdref.mat']) && po_config_data.clean_ref_ctd==1
    load([po_config_data.HISTORICAL_DIRECTORY '/historical_ctd/bad_ctdref.mat']);
else
    bad_ctd=[];
end
for ln_index = 1:m
    for ntyp = 2:4 % go through columns and check to see if this data type is supposed to be loaded
        if( ~isnan(pa_wmo_numbers(ln_index,1)) & pa_wmo_numbers(ln_index,ntyp) )
            if ntyp == 2 % the 2nd column denotes CTD data
                lo_box_data = load( strcat( po_config_data.HISTORICAL_DIRECTORY, po_config_data.HISTORICAL_CTD_PREFIX, sprintf( '%4d', pa_wmo_numbers(ln_index,1))));
                if isempty(bad_ctd)==0
                    bad_source = bad_ctd.source(bad_ctd.wmo==pa_wmo_numbers(ln_index,1));
                    if ~isempty(bad_source)
                        not_use=find(ismember(lo_box_data.source,bad_source));
                        disp('bad data ref removed')
                    else
                        not_use=[];
                    end
                else
                    bad_source=[];
                    not_use=[];
                end
                
                date_hist = changedates(lo_box_data.dates);
                lo_box_data.lat(not_use)=[];
                lo_box_data.long(not_use)=[];
                lo_box_data.dates(not_use)=[];
                
            elseif ntyp == 3 % the 3rd column denotes historical data
                lo_box_data = load( strcat( po_config_data.HISTORICAL_DIRECTORY, po_config_data.HISTORICAL_BOTTLE_PREFIX, sprintf( '%4d', pa_wmo_numbers(ln_index,1))));
            elseif ntyp == 4 % the 4th column denotes Argo data
                lo_box_data = load( strcat( po_config_data.HISTORICAL_DIRECTORY, po_config_data.HISTORICAL_ARGO_PREFIX, sprintf( '%4d', pa_wmo_numbers(ln_index,1))));
                
                % exclude Argo float being analysed from the Argo reference data selection -------
                not_use=[];
                %unique(strtok(lo_box_data.source,'_'))
                for i=1:length(lo_box_data.lat)
                    profile=lo_box_data.source{i};
                    jj=findstr(profile,'_');
                    
                    ref_float=profile(1:jj-1);
                    kk=findstr(pa_float_name, ref_float);
                    kk2=findstr('6901144', strtrim(ref_float));
                    %
                    if(isempty(kk)==0)|(isempty(kk2)==0)
                        %if(isempty(kk)==0)
                        not_use=[not_use,i];
                    end
                end
                lo_box_data.lat(not_use)=[];
                lo_box_data.long(not_use)=[];
                lo_box_data.dates(not_use)=[];
                %-----------------------------------------------------------------------
            end
            
            pa_grid_lat   = [ pa_grid_lat,   lo_box_data.lat ] ;
            pa_grid_long  = [ pa_grid_long,  lo_box_data.long ] ;
            pa_grid_dates = [ pa_grid_dates, lo_box_data.dates ] ;
        end
    end
    fclose('all');
end


% if no boxes are assigned

if( isempty( pa_grid_lat ) == 1 )
    pa_grid_lat = 999 ;
    pa_grid_long = 999 ;
    pa_grid_dates = NaN ;
end


% longitude goes from 0 to 360 degrees

ln_jj = find( pa_grid_long < 0 ) ;
pa_grid_long( ln_jj ) = 360 + pa_grid_long( ln_jj ) ;


% make sure longitude is continuous around the 0-360 degree mark

% ln_kk = find( pa_grid_long>=320 & pa_grid_long<=360 ) ;
% if( isempty( ln_kk ) == 0 )
% ln_ll = find( pa_grid_long>=0 & pa_grid_long<=40 ) ;
% pa_grid_long( ln_ll ) = 360 + pa_grid_long( ln_ll ) ;
% end
% correction ccabanes pour  pour MAIN_dmqc_speed
ln_kk1 = find( pa_grid_long>=320);
ln_kk2 =find( pa_grid_long<=40);
if( isempty( ln_kk1) == 0 & isempty(ln_kk2)==0)
    ln_ll = find( pa_grid_long>=0 & pa_grid_long<=180 ) ;
    %ln_ll = find( pa_grid_long>=0) ;
    pa_grid_long( ln_ll ) = 360 + pa_grid_long( ln_ll ) ;
end

pa_grid_dates = changedates( pa_grid_dates ) ;


% turns rows into columns

pa_grid_lat = pa_grid_lat' ;
pa_grid_long = pa_grid_long' ;


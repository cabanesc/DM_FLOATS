% -========================================================
%   USAGE :   CORR_POSITION(floatname,dacname,varargin)
%   PURPOSE : visualize and recompute and eventually correct interpolated
%   positions (QC=8)
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char) e.g.  'coriolis'
%   OPTIONNAL INPUT :
%    'METHOD' (char)  'original' (default) Dac interpolation between known position; 'sphere' interpolation 3D  adapted for high latitude
%    'CORR_NETCDF'    '1'  correction in netcdf files; '0' (default) no correction
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created C. Cabanes - 2025
% ========================================================
function CORR_POSITION(floatname,dacname,varargin)


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

% default CONFIG
PARAM.METHOD='original';
PARAM.CORR_NETCDF=0;

% Input CONFIG
if isfield(s,'METHOD')==1;PARAM.METHOD=s.METHOD;end;
if isfield(s,'CORR_NETCDF')==1;PARAM.CORR_NETCDF=s.CORR_NETCDF;end;


FLOAT_SOURCE_NETCDF=CONFIG.DIR_FTP;
%GEBCO_FILE = '/home/ccabanes/provisoire/_ressources/GEBCO_2020/GEBCO_2020.nc';
GEBCO_FILE = '/home/pharos/andro/livraisons/Livraison_0712021/soft/DecArgo_soft/_ressources/GEBCO_2020/GEBCO_2020.nc';
%GEBCO_FILE = '/export/home1/ccabanes/data/GEBCO/GEBCO_2020.nc';
for ifloat=1:length(floatname)
    thefloatname=floatname{ifloat};
    thedacname=dacname{ifloat};
    IncludeDescProf=1;
    [file_list]=select_float_files_on_ftp(thefloatname,thedacname,FLOAT_SOURCE_NETCDF,'C',IncludeDescProf);
    vertical_sampling_scheme='Primary sampling';
    Param='';
    [FLm,DimL,file_list]=create_multi_from_filelist(thefloatname,thedacname,FLOAT_SOURCE_NETCDF,file_list,vertical_sampling_scheme,Param);
    FLm=replace_fill_bynan(FLm); % add 23/01/2024
    
    M = read_netcdf_allthefile([FLOAT_SOURCE_NETCDF thedacname '/' thefloatname '/' thefloatname '_meta.nc']);
    T = read_netcdf_allthefile([FLOAT_SOURCE_NETCDF thedacname '/' thefloatname '/' thefloatname '_Rtraj.nc']);
    
    
    FLm=replace_fill_bynan(FLm);
    FLm=format_flags_char2num(FLm);
    
    % plot de la position sur le globe
    
    MAP_VISU.zone_visu=[floor(min([FLm.latitude.data;M.launch_latitude.data]))-0.2 ceil(max([FLm.latitude.data;M.launch_latitude.data]))+0.2  floor(min([FLm.longitude.data;M.launch_longitude.data]))-0.2 ceil(max([FLm.longitude.data;M.launch_longitude.data]))+0.2];
    fact=1;
    MAP_VISU.min_lat= max(floor(MAP_VISU.zone_visu(1)*fact)/fact,-90);
    MAP_VISU.max_lat= min(ceil( MAP_VISU.zone_visu(2)*fact)/fact,90);
    
    MAP_VISU.min_lon= floor(MAP_VISU.zone_visu(3)*fact)/fact;
    MAP_VISU.max_lon= ceil( MAP_VISU.zone_visu(4)*fact)/fact;
    
    % m_proj(proj,'long',[floor(zone_visu(3)) ceil(zone_visu(4))],...
    %                  'lat',[floor(zone_visu(1)) ceil(zone_visu(2))]);
    %m_grid('box','fancy','tickdir','in');
    %m_grid('tickdir','in');
    if abs(MAP_VISU.min_lat)>=70|abs(MAP_VISU.max_lat)>=80
        MAP_VISU.proj='stereographic';
        PARAM.isarctic=1;
    else
        MAP_VISU.proj='miller'; MAP_VISU.reso='LR';
        PARAM.isarctic=0;
    end
    %reso='LR';proj='miller';
    
    %     [hf,ha]=fct_pltmap_traj(MAP_VISU,PARAM,1000,2000);
    %
    %     m_scatter(FLm.longitude.data,FLm.latitude.data,40,'^w','filled');
    %     hold on
    %    % hj1=m_scatter(M.launch_longitude.data,M.launch_latitude.data,60,'sr','filled');
    %
    %     isPosQc1=find(FLm.position_qc.data==1);
    %     hj2=m_scatter(FLm.longitude.data(isPosQc1),FLm.latitude.data(isPosQc1),60,'^g','filled');
    %     m_plot(FLm.longitude.data(isPosQc1),FLm.latitude.data(isPosQc1))
    %     hj1=m_scatter(M.launch_longitude.data,M.launch_latitude.data,60,'sr','filled');
    %     load ('bathy_28_colormap.mat');
    %     colormap(newmap4)
    %     colorbar
    % Interpole entre deux positions qc 1, aux dates juld
    
    % prépare les données
    floatDataAll.cycleNumber=FLm.cycle_number.data';
    floatDataAll.direction=FLm.direction.data';
    floatDataAll.latitude=FLm.latitude.data';
    floatDataAll.longitude=FLm.longitude.data';
    floatDataAll.time=FLm.juld_location.data';
    floatDataAll.profPresMax=max(FLm.pres.data(:,:)');
    floatDataAll.positionQc=FLm.position_qc.data(:)';
    
    T=replace_fill_bynan(T);
    
    floatDataAll.grounded=NaN*ones(length(floatDataAll.cycleNumber),1);
    
    for idCy = 1:length(FLm.cycle_number.data)
        cyNum = FLm.cycle_number.data(idCy);
        idForCy = find(floatDataAll.cycleNumber == cyNum);
        groundedval=T.grounded.data(T.cycle_number_index.data == cyNum);
        
        if ~(isempty(groundedval))&~(isempty(idForCy))
            if groundedval=='Y'
                floatDataAll.grounded(idForCy)=1;
            elseif groundedval=='N'
                floatDataAll.grounded(idForCy)=0;
            end
        end
    end
    
    
    % retrieve GEBCO depth
    floatDataAll.gebcoDepth = get_gebco_depth(floatDataAll.longitude, floatDataAll.latitude, GEBCO_FILE);
    floatDataAll.index=[1:length(floatDataAll.cycleNumber)];
    
    floatDataAll.time(isnan(floatDataAll.time))=FLm.juld.data(isnan(floatDataAll.time))';
    % garde seulement les points avec position_qc de 1 2, 9(mais date connue) ou 8
    isKeepQC = (floatDataAll.positionQc==1)|(floatDataAll.positionQc==2)|(floatDataAll.positionQc==8)|(floatDataAll.positionQc==9&~isnan(floatDataAll.time));
    floatDataAll.cycleNumber = floatDataAll.cycleNumber(isKeepQC);
    floatDataAll.direction = floatDataAll.direction(isKeepQC);
    floatDataAll.latitude = floatDataAll.latitude(isKeepQC);
    floatDataAll.longitude = floatDataAll.longitude(isKeepQC);
    floatDataAll.time = floatDataAll.time(isKeepQC);
    floatDataAll.profPresMax = floatDataAll.profPresMax(isKeepQC);
    floatDataAll.positionQc = floatDataAll.positionQc(isKeepQC);
    floatDataAll.grounded = floatDataAll.grounded(isKeepQC);
    floatDataAll.gebcoDepth =  floatDataAll.gebcoDepth(isKeepQC);
    floatDataAll.index = floatDataAll.index(isKeepQC);
    
    % stock des résultats
    floatDataAll.latInterp = NaN*floatDataAll.latitude;
    floatDataAll.lonInterp = NaN*floatDataAll.longitude;
    
    segments = find_missing_segments(floatDataAll);
    
    % Traiter chaque segment manquant
    for i = 1:size(segments, 1)
        start_idx = segments(i, 1);
        end_idx = segments(i, 2);
        
        fprintf('   - Traitement segment %d: cycles %d à %d\n', i, floatDataAll.cycleNumber(start_idx), floatDataAll.cycleNumber(end_idx));
        
        
        
        % Calculer les indices à extraire du résultat
        segment_data = extract_segment_data(floatDataAll,M, start_idx, end_idx);
        missing_start_rel = segment_data.missing_start_rel;
        missing_end_rel = segment_data.missing_end_rel;
        
        juld_query = segment_data.time(missing_start_rel:missing_end_rel);
        % Interpolation lineaire sur la sphere
        if strcmp(PARAM.METHOD,'sphere')
            [lat_interp, lon_interp] = interpolate_geopos_sphere_improved(segment_data.latitude(1), segment_data.longitude(1),segment_data.latitude(end), segment_data.longitude(end), juld_query, segment_data.time(1), segment_data.time(end));
        elseif strcmp(PARAM.METHOD,'original')
            [lat_interp, lon_interp] = interpolate_between_2_locations(segment_data.time(1),segment_data.longitude(1), segment_data.latitude(1),segment_data.time(end), segment_data.longitude(end), segment_data.latitude(end), juld_query);
        elseif strcmp(PARAM.METHOD,'bathy')
            [lat_interp, lon_interp, depth_path] = interpolate_along_bathymetry(segment_data.latitude(1), segment_data.longitude(1),segment_data.latitude(end), segment_data.longitude(end), juld_query, segment_data.time(1), segment_data.time(end),'n_waypoints',length(juld_query));
        else
            error ('interpolation method unknown')
        end
        
        %         figure(hf)
        %         hj3=m_scatter(lon_interp,lat_interp,40,'^k','filled');
        %         legend([hj1;hj2; hj3],{'Launch position', 'POS_QC=1','interpolated 3D linear'},'interpreter','none')
        
        floatDataAll.latInterp(segment_data.target_start:segment_data.target_end)=lat_interp;
        floatDataAll.lonInterp(segment_data.target_start:segment_data.target_end)=lon_interp;
        
        %         % lecture des fichiers d'estimation ( terrain following)
        %         file_estimate=[CONFIG.DIR_PLOT 'under_ice/' thefloatname '/estimate_profile_locations_' thefloatname '.csv'];
        %
        %         if exist(file_estimate)
        %             results = readtable(file_estimate, 'NumHeaderLines', 0);
        %             hj4=m_scatter(results.TrajLon,results.TrajLat,'*r');
        %            % legend([hj1;hj2; hj3;hj4],{'Launch position', 'POS_QC=1','interpolated 3D linear','estimation terrain_following'},'interpreter','none')
        %         end
        %
        %         file_estimate=[CONFIG.DIR_PLOT 'under_ice/' thefloatname '/estimate_profile_locations_' thefloatname '_highLat.csv'];
        %
        %         if exist(file_estimate)
        %             results = readtable(file_estimate, 'NumHeaderLines', 0);
        %             hj5=m_scatter(results.TrajLon,results.TrajLat,'om');
        %           %  legend([hj1;hj2; hj3;hj4; hj5],{'Launch position', 'POS_QC=1','interpolated 3D linear','estimation terrain_following','estimation terrain_following High Lat'},'interpreter','none')
        %         end
        
        
        
        
    end
    close all
    
    % tracé des trajectoires
    [hf,ha]=fct_pltmap_traj(MAP_VISU,PARAM,1000,2000);
    
    m_scatter(FLm.longitude.data,FLm.latitude.data,40,'^w','filled');
    hold on
    % hj1=m_scatter(M.launch_longitude.data,M.launch_latitude.data,60,'sr','filled');
    
    isPosQc1=find(FLm.position_qc.data==1);
    hj2=m_scatter(FLm.longitude.data(isPosQc1),FLm.latitude.data(isPosQc1),60,'^g','filled');
    m_plot(FLm.longitude.data(isPosQc1),FLm.latitude.data(isPosQc1))
    hj1=m_scatter(M.launch_longitude.data,M.launch_latitude.data,60,'sr','filled');
    load ('bathy_28_colormap.mat');
    colormap(newmap4)
    colorbar
    hj3=m_scatter(floatDataAll.lonInterp,floatDataAll.latInterp,40,'^k','filled');
    %m_scatter(floatDataAll.lonInterp(floatDataAll.grounded==1),floatDataAll.latInterp(floatDataAll.grounded==1),40,'^m');
    legend([hj1;hj2; hj3],{'Launch position', 'POS_QC=1','interpolated 3D linear'},'interpreter','none')
    % lecture des fichiers d'estimation ( terrain following)
    file_estimate=[CONFIG.DIR_PLOT 'under_ice/' thefloatname '/estimate_profile_locations_' thefloatname '.csv'];
    
    if exist(file_estimate)
        results = readtable(file_estimate, 'NumHeaderLines', 0);
        hj4=m_scatter(results.TrajLon,results.TrajLat,'*r');
        % legend([hj1;hj2; hj3;hj4],{'Launch position', 'POS_QC=1','interpolated 3D linear','estimation terrain_following'},'interpreter','none')
    end
    
    % test de la validité des interpolations
    if sum(floatDataAll.positionQc==8|floatDataAll.positionQc==9)~=sum(~isnan(floatDataAll.latInterp))|sum(floatDataAll.positionQc==8|floatDataAll.positionQc==9)~=sum(~isnan(floatDataAll.lonInterp))
        error('Toutes les données avec positions manquantes n''ont pas pu etre interpolees')
    end
    
    figure(2)
    floatDataAll.gebcoDepthNew = get_gebco_depth(floatDataAll.lonInterp,floatDataAll.latInterp, GEBCO_FILE);
    %
    h1=plot(floatDataAll.cycleNumber,floatDataAll.profPresMax,'+');
    hold on
    h2=plot(floatDataAll.cycleNumber,floatDataAll.gebcoDepthNew,'+');
    set(gca,'Ydir','reverse')
    h3=plot(floatDataAll.cycleNumber(floatDataAll.grounded==1),floatDataAll.profPresMax(floatDataAll.grounded==1),'om');
    grid on
    box on
    xlabel('cycle Number')
    ylabel('pressure(db)')
    legend([h1;h2; h3],{'Maximal Pressure of the profile', 'Bathymetry','grounded'},'interpreter','none');
    
    if PARAM.CORR_NETCDF==1
    % correction des positions dans les fichiers netcdf
        for ifiles=1:length(file_list)
    
            file_name = [CONFIG.DIR_FTP thedacname '/' thefloatname '/profiles/' file_list{ifiles} ];
            [F,Dim,G]=read_netcdf_allthefile(file_name);
            F=replace_fill_bynan(F);
            isposQC8 = findstr_tab(F.position_qc.data,'8')|findstr_tab(F.position_qc.data,'9');
            ismodif=0;
            for iprof=1:length(F.cycle_number.data)
                if isposQC8(iprof)
                    % remplace la position par la nouvelle interpolation  floatDataAll
                    isInd = (floatDataAll.cycleNumber==F.cycle_number.data(iprof))&(floatDataAll.direction==F.direction.data(iprof))&~isnan(floatDataAll.latInterp)&~isnan(floatDataAll.lonInterp);
                    %verifie que la nouvelle poisition n'est pas sur terre
                    gebcoDepth = get_gebco_depth(floatDataAll.lonInterp(isInd), floatDataAll.latInterp(isInd), GEBCO_FILE);
    
                    if sum(isInd)==1&gebcoDepth>0
                        F.latitude.data(iprof)=floatDataAll.latInterp(isInd);
                        F.longitude.data(iprof)=floatDataAll.lonInterp(isInd);
                        F.position_qc.data(iprof)='8';
                        ismodif=1;
                        disp(['Cycle: ' num2str(F.cycle_number.data(iprof)) ': position interpolée'])
                    else
                       F.latitude.data(iprof)=NaN;
                       F.longitude.data(iprof)=NaN;
                       F.position_qc.data(iprof)='9';
                       disp(['Cycle: ' num2str(F.cycle_number.data(iprof)) ': position manquante'])
                    end
                end
            end
            % sauvegarde des fichiers corrigés
            if ismodif==1
                F=replace_nan_byfill(F);
                create_netcdf_allthefile(F,Dim,file_name,G)
            end
            % F.latitude.data
        end
    end
    
    
end

end

%------------------------------------------------------------------------
%  FONCTIONS NECESSAIRES
% ------------------------------------------------------------------------

function [lat_interp, lon_interp] = interpolate_geopos_sphere(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end)
%INTERPOLATE_GEOPOS_SPHERE Interpolation géodésique entre deux points
% sur une sphère sans toolbox, basée sur l'interpolation linéaire en 3D.
%
% Entrées :
%   lat1, lon1 : Point de départ (en degrés)
%   lat2, lon2 : Point d'arrivée (en degrés)
%   juld_query : Dates (jours julien) où on veut les positions
%   juld_start : Début du trajet
%   juld_end   : Fin du trajet
%
% Sorties :
%   lat_interp, lon_interp : Positions interpolées (en degrés)

% Rayon terrestre (approx.) en km
R = 6371;

% Conversion en radians
lat1 = deg2rad(lat1); lon1 = deg2rad(lon1);
lat2 = deg2rad(lat2); lon2 = deg2rad(lon2);

% Vecteurs 3D sur la sphère
p1 = R * [cos(lat1)*cos(lon1), cos(lat1)*sin(lon1), sin(lat1)];
p2 = R * [cos(lat2)*cos(lon2), cos(lat2)*sin(lon2), sin(lat2)];

% Angle entre les deux points (en radians)
omega = acos(dot(p1, p2) / (R^2));

% Interpolation fractionnelle
f = (juld_query - juld_start) / (juld_end - juld_start);

lat_interp = zeros(size(f));
lon_interp = zeros(size(f));

for i = 1:length(f)
    if omega == 0
        p = p1;
    else
        % Interpolation de type "slerp" (spherical linear interpolation)
        A = sin((1 - f(i)) * omega) / sin(omega);
        B = sin(f(i) * omega) / sin(omega);
        p = A * p1 + B * p2;
    end
    
    % Conversion inverse en lat/lon
    x = p(1); y = p(2); z = p(3);
    lat_interp(i) = rad2deg(asin(z / R));
    lon_interp(i) = rad2deg(atan2(y, x));
end

end


function [lat_interp, lon_interp] = interpolate_geopos(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end)
%INTERPOLATE_GEOPOS Interpolation lineaire entre deux points (lon,lat)
%
% Entrées :
%   lat1, lon1 : Point de départ (en degrés)
%   lat2, lon2 : Point d'arrivée (en degrés)
%   juld_query : Dates (jours julien) où on veut les positions
%   juld_start : Début du trajet
%   juld_end   : Fin du trajet
%
% Sorties :
%   lat_interp, lon_interp : Positions interpolées (en degrés)

% Rayon terrestre (approx.) en km
R = 6371;

% Conversion en radians
lat1 = deg2rad(lat1); lon1 = deg2rad(lon1);
lat2 = deg2rad(lat2); lon2 = deg2rad(lon2);

% Vecteurs 3D sur la sphère
p1 = R * [cos(lat1)*cos(lon1), cos(lat1)*sin(lon1), sin(lat1)];
p2 = R * [cos(lat2)*cos(lon2), cos(lat2)*sin(lon2), sin(lat2)];

% Angle entre les deux points (en radians)
omega = acos(dot(p1, p2) / (R^2));

% Interpolation fractionnelle
f = (juld_query - juld_start) / (juld_end - juld_start);

lat_interp = zeros(size(f));
lon_interp = zeros(size(f));

for i = 1:length(f)
    if omega == 0
        p = p1;
    else
        % Interpolation de type "slerp" (spherical linear interpolation)
        A = sin((1 - f(i)) * omega) / sin(omega);
        B = sin(f(i) * omega) / sin(omega);
        p = A * p1 + B * p2;
    end
    
    % Conversion inverse en lat/lon
    x = p(1); y = p(2); z = p(3);
    lat_interp(i) = rad2deg(asin(z / R));
    lon_interp(i) = rad2deg(atan2(y, x));
end

end
% ------------------------------------------------------------------------------
% Get GEBCO depth associated to a list of locations.
%
% SYNTAX :
%  [o_depth] = get_gebco_depth(a_lon, a_lat, a_gebcoPathFileName)
%
% INPUT PARAMETERS :
%   a_lon               : latitudes
%   a_lat               : longitudes
%   a_gebcoPathFileName : GEBCO path file name
%
% OUTPUT PARAMETERS :
%   o_depth : depths
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   03/01/2022 - RNU - creation
% ------------------------------------------------------------------------------
function [o_depth] = get_gebco_depth(a_lon, a_lat, a_gebcoPathFileName)

[elevOri] = get_gebco_elev_point(a_lon, a_lat, a_gebcoPathFileName);
elev = mean(elevOri, 2);
if (any(isnan(elev)))
    idNan = find(isnan(elev));
    for idL = idNan'
        elev(idL) = mean(elevOri(idL, ~isnan(elevOri(idL, :))));
    end
end
o_depth = -elev';

end
% ------------------------------------------------------------------------------
% Retrieve the surounding elevations of a list of locations from the GEBCO 2019
% file.
%
% SYNTAX :
%  [o_elev] = get_gebco_elev_point(a_lon, a_lat, a_gebcoFileName)
%
% INPUT PARAMETERS :
%   a_lon           : list of location longitudes
%   a_lat           : list of location atitudes
%   a_gebcoFileName : GEBCO 2019 file path name
%
% OUTPUT PARAMETERS :
%   o_elev : surrounding elevations of each location
%            (size(o_elev) = [length(a_lon) 4]
%             4 elevations are generally provided [elevSW elevNW elevSE elevNE]
%             when only 1 or 2 are provided other ones are set to NaN)
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   04/29/2020 - RNU - creation
% ------------------------------------------------------------------------------
function [o_elev] = get_gebco_elev_point(a_lon, a_lat, a_gebcoFileName)

% output parameters initialization
o_elev = nan(length(a_lon), 4);


% check inputs
if (a_lon < -180)
    fprintf('ERROR: get_gebco_elev_point: input lon < -180\n');
    return
end
if (a_lon >= 360)
    fprintf('ERROR: get_gebco_elev_point: input lon >= 360\n');
    return
end
if (a_lat < -90)
    fprintf('ERROR: get_gebco_elev_point: input lat < -90\n');
    return
elseif (a_lat > 90)
    fprintf('ERROR: get_gebco_elev_point: input lat > 90\n');
    return
end

if (a_lon >= 180)
    a_lon = a_lon - 360;
end

% check GEBCO file exists
if ~(exist(a_gebcoFileName, 'file') == 2)
    fprintf('ERROR: GEBCO file not found (%s)\n', a_gebcoFileName);
    return
end

% open NetCDF file
fCdf = netcdf.open(a_gebcoFileName, 'NC_NOWRITE');
if (isempty(fCdf))
    fprintf('RTQC_ERROR: Unable to open NetCDF input file: %s\n', a_gebcoFileName);
    return
end

lonVarId = netcdf.inqVarID(fCdf, 'lon');
latVarId = netcdf.inqVarID(fCdf, 'lat');
elevVarId = netcdf.inqVarID(fCdf, 'elevation');

lon = netcdf.getVar(fCdf, lonVarId);
lat = netcdf.getVar(fCdf, latVarId);
minLon = min(lon);
maxLon = max(lon);

for idP = 1:length(a_lat)
    
    if (isnan(a_lat(idP)) || isnan(a_lon(idP)))
        continue
    end
    
    idLigStart = find(lat <= a_lat(idP), 1, 'last');
    if (isempty(idLigStart))
        idLigStart = 1;
    end
    idLigEnd = find(lat >= a_lat(idP), 1, 'first');
    if (isempty(idLigEnd))
        idLigEnd = length(lat);
    end
    %    latVal = lat(fliplr(idLigStart:idLigEnd));
    
    % a_lon(idP) is in the [-180, 180[ interval
    % it can be in 3 zones:
    % case 1: [-180, minLon[
    % case 2: [minLon, maxLon]
    % case 3: ]maxLon, -180[
    if ((a_lon(idP) >= minLon) && (a_lon(idP) <= maxLon))
        % case 2
        idColStart = find(lon <= a_lon(idP), 1, 'last');
        idColEnd = find(lon >= a_lon(idP), 1, 'first');
        
        elev = nan(length(idLigStart:idLigEnd), length(idColStart:idColEnd));
        for idL = idLigStart:idLigEnd
            elev(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 idColStart-1]), fliplr([1 length(idColStart:idColEnd)]))';
        end
        
        %       lonVal = lon(idColStart:idColEnd);
    elseif (a_lon(idP) < minLon)
        % case 1
        elev1 = nan(length(idLigStart:idLigEnd), 1);
        for idL = idLigStart:idLigEnd
            elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 length(lon)-1]), fliplr([1 1]))';
        end
        
        %       lonVal1 = lon(end);
        
        elev2 = nan(length(idLigStart:idLigEnd), 1);
        for idL = idLigStart:idLigEnd
            elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 0]), fliplr([1 1]))';
        end
        
        %       lonVal2 = lon(1) + 360;
        
        elev = cat(2, elev1, elev2);
        %       lonVal = cat(1, lonVal1, lonVal2);
        clear elev1 elev2
    elseif (a_lon(idP) > maxLon)
        % case 3
        elev1 = nan(length(idLigStart:idLigEnd), 1);
        for idL = idLigStart:idLigEnd
            elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 length(lon)-1]), fliplr([1 1]))';
        end
        
        %       lonVal1 = lon(end);
        
        elev2 = nan(length(idLigStart:idLigEnd), 1);
        for idL = idLigStart:idLigEnd
            elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 0]), fliplr([1 1]))';
        end
        
        %       lonVal2 = lon(1) + 360;
        
        elev = cat(2, elev1, elev2);
        %       lonVal = cat(1, lonVal1, lonVal2);
        clear elev1 elev2
    end
    
    if (~isempty(elev))
        if (size(elev, 1) == 2)
            if (size(elev, 2) == 2)
                o_elev(idP, 1) = elev(2, 1);
                o_elev(idP, 2) = elev(1, 1);
                o_elev(idP, 3) = elev(2, 2);
                o_elev(idP, 4) = elev(1, 2);
            else
                o_elev(idP, 1) = elev(2);
                o_elev(idP, 2) = elev(1);
            end
        else
            if (size(elev, 2) == 2)
                o_elev(idP, 1) = elev(1, 1);
                o_elev(idP, 3) = elev(1, 2);
            else
                o_elev(idP, 1) = elev;
            end
        end
    end
    
    clear elev
end

netcdf.close(fCdf);

clear lon lat

end




function segments = find_missing_segments(floatData)
% Trouve les segments continus de positions manquantes

missing_indices = find(floatData.positionQc == 8|floatData.positionQc == 9);

if isempty(missing_indices)
    segments = [];
    return;
end

% Grouper les indices consécutifs
segments = [];
start_idx = missing_indices(1);

for i = 2:length(missing_indices)
    if missing_indices(i) ~= missing_indices(i-1) + 1
        % Fin d'un segment
        segments = [segments; start_idx, missing_indices(i-1)];
        start_idx = missing_indices(i);
    end
end

% Dernier segment
segments = [segments; start_idx, missing_indices(end)];
end



function segment_data = extract_segment_data(floatData,M, start_idx, end_idx)
% Extrait les données pour le segment (avec points de départ/arrivée)

% Pour l'estimation, on a besoin du point avant le segment manquant et du point après
% Donc on inclut start_idx-1 et end_idx+1 comme points connus
actual_start = max(1, start_idx - 1);
actual_end = min(length(floatData.latitude), end_idx + 1);
segment_data = struct();
while (actual_start~=0)&&(floatData.positionQc(actual_start)~=1&&floatData.positionQc(actual_start)~=2)
    actual_start = actual_start - 1;
end
while (actual_end~=length(floatData.latitude))&&(floatData.positionQc(actual_end)~=1&&floatData.positionQc(actual_end)~=2)
    actual_end = actual_end + 1;
    % (actual_end~=length(floatData.latitude))
    % (floatData.positionQc(actual_end)~=1||floatData.positionQc(actual_end)~=2)
    % floatData.positionQc(actual_end)
end
if actual_start==0
    disp('Utilisation de la position/date de mise à l''eau')
    actual_start=1;
    %start_idx=start_idx+1;
    %end_idx=min(length(floatData.latitude), end_idx + 1);
    juldLaunch=datenum(M.launch_date.data','yyyymmddHHMMSS')-datenum('01011950','ddmmyyyy');
    segment_data.cycleNumber =[0,floatData.cycleNumber(actual_start:actual_end)];
    segment_data.latitude = [M.launch_latitude.data,floatData.latitude(actual_start:actual_end)];
    segment_data.longitude = [M.launch_longitude.data,floatData.longitude(actual_start:actual_end)];
    segment_data.time = [juldLaunch,floatData.time(actual_start:actual_end)];
    segment_data.profPresMax = [NaN,floatData.profPresMax(actual_start:actual_end)];
    segment_data.grounded = [NaN,floatData.grounded(actual_start:actual_end)'];
    segment_data.gebcoDepth = [NaN,floatData.gebcoDepth(actual_start:actual_end)];
    
    segment_data.target_start = start_idx;  % Premier indice des positions à estimer
    segment_data.target_end = end_idx;      % Dernier indice des positions à estimer
    % % Calculer les indices relatifs dans le segment
    segment_data.missing_start_rel = 2;  % Position dans segment_data
    segment_data.missing_end_rel = end_idx -start_idx + 2;      % Position dans segment_data
else
    segment_data.cycleNumber =floatData.cycleNumber(actual_start:actual_end);
    segment_data.latitude = floatData.latitude(actual_start:actual_end);
    segment_data.longitude = floatData.longitude(actual_start:actual_end);
    segment_data.time = floatData.time(actual_start:actual_end);
    segment_data.profPresMax = floatData.profPresMax(actual_start:actual_end);
    segment_data.grounded = floatData.grounded(actual_start:actual_end);
    segment_data.gebcoDepth = floatData.gebcoDepth(actual_start:actual_end);
    
    segment_data.target_start = start_idx;  % Premier indice des positions à estimer
    segment_data.target_end = end_idx;      % Dernier indice des positions à estimer
     % % Calculer les indices relatifs dans le segment
    segment_data.missing_start_rel = start_idx - actual_start + 1;  % Position dans segment_data
    segment_data.missing_end_rel = end_idx - actual_start + 1;      % Position dans segment_data
    
end

% Stocker les indices pour le mapping
% segment_data.actual_start = actual_start;
% segment_data.actual_end = actual_end;



end

function data = read_profile_csv(filename)
%READ_PROFILE_CSV Lit correctement le fichier CSV avec détection des types
% et conversion automatique des cellules

% Lire la table avec détection d'en-tête
opts = detectImportOptions(filename, 'Delimiter', ';', 'NumHeaderLines', 0);
opts = setvaropts(opts, opts.VariableNames, 'WhitespaceRule', 'preserve');
opts = setvaropts(opts, opts.VariableNames, 'EmptyFieldRule', 'auto');

% Forcer la 1re ligne comme noms de colonnes
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];  % Commence à lire les données à la ligne 2

% Lire la table
T = readtable(filename, opts);

% Corriger les colonnes texte contenant des dates
for i = 1:width(T)
    if iscellstr(T{:, i}) || isstring(T{:, i})
        try
            T{:, i} = datetime(T{:, i}, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
        catch
            % Laisser tel quel si la conversion échoue
        end
    end
end

% Convertir en structure pour accès facile
data = table2struct(T, 'ToScalar', true);
end


function [lat_interp, lon_interp] = interpolate_geopos_sphere_improved(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end)
%INTERPOLATE_GEOPOS_SPHERE_IMPROVED Interpolation géodésique entre deux points
% sur une sphère sans toolbox, basée sur l'interpolation SLERP optimisée.
%
% Cette version corrige plusieurs problèmes numériques et ajoute des
% optimisations pour les cas particuliers.
%
% Entrées :
% lat1, lon1 : Point de départ (en degrés)
% lat2, lon2 : Point d'arrivée (en degrés)
% juld_query : Dates (jours julien) où on veut les positions
% juld_start : Début du trajet
% juld_end : Fin du trajet
%
% Sorties :
% lat_interp, lon_interp : Positions interpolées (en degrés)

%% Validation des entrées
if juld_end <= juld_start
    error('juld_end doit être supérieur à juld_start');
end

% Gérer les cas où juld_query est en dehors de la plage
juld_query = max(juld_start, min(juld_end, juld_query));

%% Constantes
R = 6371; % Rayon terrestre (approx.) en km

%% Conversion en radians
lat1_rad = deg2rad(lat1);
lon1_rad = deg2rad(lon1);
lat2_rad = deg2rad(lat2);
lon2_rad = deg2rad(lon2);

%% Gestion des longitudes pour éviter les problèmes d'antimeridien
% Si la différence de longitude > 180°, utiliser le chemin le plus court
dlon = lon2_rad - lon1_rad;
if abs(dlon) > pi
    if dlon > 0
        lon1_rad = lon1_rad + 2*pi;
    else
        lon2_rad = lon2_rad + 2*pi;
    end
end

%% Vecteurs 3D sur la sphère unitaire (normalisés)
p1 = [cos(lat1_rad)*cos(lon1_rad), cos(lat1_rad)*sin(lon1_rad), sin(lat1_rad)];
p2 = [cos(lat2_rad)*cos(lon2_rad), cos(lat2_rad)*sin(lon2_rad), sin(lat2_rad)];

%% Calcul de l'angle entre les deux points
dot_product = dot(p1, p2);
% Clamp pour éviter les erreurs numériques avec acos
dot_product = max(-1, min(1, dot_product));
omega = acos(dot_product);

%% Interpolation fractionnelle
f = (juld_query - juld_start) / (juld_end - juld_start);

%% Pré-allocation pour optimiser
lat_interp = zeros(size(f));
lon_interp = zeros(size(f));

%% Cas particuliers pour optimiser
tolerance = 1e-10;

if omega < tolerance
    % Points très proches ou identiques : interpolation linéaire simple
    lat_interp = lat1 + f * (lat2 - lat1);
    lon_interp = lon1 + f * (lon2 - lon1);
    
    % Normaliser les longitudes dans [-180, 180]
    lon_interp = mod(lon_interp + 180, 360) - 180;
    
elseif abs(omega - pi) < tolerance
    % Points antipodaux : trajectoire ambiguë, utiliser interpolation linéaire
    warning('Points antipodaux détectés, interpolation linéaire utilisée');
    lat_interp = lat1 + f * (lat2 - lat1);
    lon_interp = lon1 + f * (lon2 - lon1);
    
    % Normaliser les longitudes
    lon_interp = mod(lon_interp + 180, 360) - 180;
    
else
    % Interpolation SLERP normale
    sin_omega = sin(omega);
    
    for i = 1:length(f)
        if f(i) <= 0
            % Avant le début
            p = p1;
        elseif f(i) >= 1
            % Après la fin
            p = p2;
        else
            % SLERP classique
            A = sin((1 - f(i)) * omega) / sin_omega;
            B = sin(f(i) * omega) / sin_omega;
            p = A * p1 + B * p2;
            
            % Re-normaliser pour maintenir la précision numérique
            p = p / norm(p);
        end
        
        % Conversion inverse en lat/lon
        x = p(1); y = p(2); z = p(3);
        
        % Clamp z pour éviter les erreurs numériques avec asin
        z = max(-1, min(1, z));
        
        lat_interp(i) = rad2deg(asin(z));
        lon_interp(i) = rad2deg(atan2(y, x));
    end
end

%% Normalisation finale des longitudes dans [-180, 180]
lon_interp = mod(lon_interp + 180, 360) - 180;

end





%% ========================================================================
%% INTERPOLATION LINÉAIRE DE BASE (SUR LE PLAN)
%% ========================================================================
function [lat_interp, lon_interp] = interpolate_geopos_linear_plane(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end)
%INTERPOLATE_GEOPOS_LINEAR_PLANE Interpolation linéaire simple sur le plan
%
% Interpolation linéaire basique en traitant longitude et latitude comme
% des coordonnées cartésiennes.
%
% Entrées :
% lat1, lon1 : Point de départ (en degrés)
% lat2, lon2 : Point d'arrivée (en degrés)
% juld_query : Dates (jours julien) où on veut les positions
% juld_start : Début du trajet
% juld_end : Fin du trajet
%
% Sorties :
% lat_interp, lon_interp : Positions interpolées (en degrés)

%% Validation des entrées
if juld_end <= juld_start
    error('juld_end doit être supérieur à juld_start');
end

%% Facteur d'interpolation
f = (juld_query - juld_start) / (juld_end - juld_start);

%% Gestion de l'antimeridien pour éviter les trajectoires aberrantes
dlon = lon2 - lon1;
if abs(dlon) > 180
    if dlon > 0
        lon1 = lon1 + 360;  % Décaler lon1
    else
        lon2 = lon2 + 360;  % Décaler lon2
    end
end

%% Interpolation linéaire directe
lat_interp = lat1 + f .* (lat2 - lat1);
lon_interp = lon1 + f .* (lon2 - lon1);

%% Normalisation des longitudes dans [-180, 180]
lon_interp = mod(lon_interp + 180, 360) - 180;

%% Clamp les latitudes dans [-90, 90] (sécurité)
lat_interp = max(-90, min(90, lat_interp));

end


function [o_interpLocLat, o_interpLocLon] = interpolate_between_2_locations(...
    a_firstLocDate, a_firstLocLon, a_firstLocLat, ...
    a_secondLocDate, a_secondLocLon, a_secondLocLat, ...
    a_interpDate)

% output parameters initialization
o_interpLocLon = [];
o_interpLocLat = [];


% interpolate between the locations
if (((abs(a_firstLocLon) > 90) && (abs(a_secondLocLon) > 90)) && ...
        (((a_firstLocLon > 0) && (a_secondLocLon < 0)) || ((a_secondLocLon > 0) && (a_firstLocLon < 0))))
    % the float crossed the date line
    if (a_secondLocLon < 0)
        a_secondLocLon = a_secondLocLon + 360;
    else
        a_firstLocLon = a_firstLocLon + 360;
    end
    o_interpLocLon = interp1q([a_firstLocDate; a_secondLocDate], [a_firstLocLon; a_secondLocLon], a_interpDate');
    if (o_interpLocLon >= 180)
        o_interpLocLon = o_interpLocLon - 360;
    end
else
    o_interpLocLon = interp1q([a_firstLocDate; a_secondLocDate], [a_firstLocLon; a_secondLocLon], a_interpDate');
end
o_interpLocLat = interp1q([a_firstLocDate; a_secondLocDate], [a_firstLocLat; a_secondLocLat], a_interpDate');

end

%% ========================================================================
%% INTERPOLATION QUI SUIT LES ISOBATHS (SIMPLE)
%% ========================================================================

function [lat_interp, lon_interp, depth_path] = interpolate_along_bathymetry(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end, varargin)
%INTERPOLATE_ALONG_BATHYMETRY Interpolation suivant la bathymétrie avec m_tbase
%
% Cette fonction interpole entre deux positions en suivant les contours
% bathymétriques plutôt qu'une ligne droite, ce qui est plus réaliste
% pour les flotteurs Argo qui dérivent à profondeur constante.
%
% SYNTAXE :
%   [lat_interp, lon_interp] = interpolate_along_bathymetry(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end)
%   [lat_interp, lon_interp, depth_path] = interpolate_along_bathymetry(..., 'param', value, ...)
%
% ENTRÉES :
%   lat1, lon1   : Position de départ (degrés)
%   lat2, lon2   : Position d'arrivée (degrés)
%   juld_query   : Dates où calculer les positions
%   juld_start   : Date de départ
%   juld_end     : Date d'arrivée
%
% PARAMÈTRES OPTIONNELS :
%   'method'           : 'isobath', 'gradient', 'constrained' (défaut: 'isobath')
%   'target_depth'     : Profondeur cible en mètres (défaut: profondeur moyenne)
%   'depth_tolerance'  : Tolérance sur la profondeur en mètres (défaut: 200m)
%   'n_waypoints'      : Nombre de points de passage (défaut: 20)
%   'max_iterations'   : Nb max d'itérations pour optimisation (défaut: 100)
%   'fallback'         : Méthode de repli si échec ('linear', 'sphere', défaut: 'sphere')
%
% SORTIES :
%   lat_interp   : Latitudes interpolées
%   lon_interp   : Longitudes interpolées
%   depth_path   : Profondeurs le long du chemin (optionnel)

%% Validation des entrées et paramètres par défaut
p = inputParser;
addRequired(p, 'lat1', @isnumeric);
addRequired(p, 'lon1', @isnumeric);
addRequired(p, 'lat2', @isnumeric);
addRequired(p, 'lon2', @isnumeric);
addRequired(p, 'juld_query', @isnumeric);
addRequired(p, 'juld_start', @isnumeric);
addRequired(p, 'juld_end', @isnumeric);

addParameter(p, 'method', 'isobath', @(x) ismember(x, {'isobath', 'gradient', 'constrained'}));
addParameter(p, 'target_depth', [], @isnumeric);
addParameter(p, 'depth_tolerance', 200, @isnumeric);
addParameter(p, 'n_waypoints', 20, @isnumeric);
addParameter(p, 'max_iterations', 100, @isnumeric);
addParameter(p, 'fallback', 'sphere', @(x) ismember(x, {'linear', 'sphere'}));

parse(p, lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end, varargin{:});
params = p.Results;

%% Vérification que m_tbase est disponible
if ~exist('m_tbase', 'file')
    warning('m_tbase non disponible. Utilisation de la méthode de repli: %s', params.fallback);
    if strcmp(params.fallback, 'linear')
        [lat_interp, lon_interp] = interpolate_linear_fallback(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
    else
        [lat_interp, lon_interp] = interpolate_sphere_fallback(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
    end
    depth_path = [];
    return;
end

%% Calcul de la grille bathymétrique
try
    [lat_grid, lon_grid, depth_grid] = get_bathymetry_grid(lat1, lon1, lat2, lon2);
    
    if isempty(depth_grid)
        warning('Impossible d''obtenir la bathymétrie. Méthode de repli utilisée.');
        [lat_interp, lon_interp] = fallback_interpolation(params.fallback, lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
        depth_path = [];
        return;
    end
    
catch ME
    warning('Erreur avec m_tbase: %s. Méthode de repli utilisée.');
    [lat_interp, lon_interp] = fallback_interpolation(params.fallback, lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
    depth_path = [];
    return;
end

%% Détermination de la profondeur cible
if isempty(params.target_depth)
    % Profondeur cible = moyenne des profondeurs aux points de départ et d'arrivée
    depth1 = interp2(lon_grid, lat_grid, depth_grid, lon1, lat1, 'linear', NaN);
    depth2 = interp2(lon_grid, lat_grid, depth_grid, lon2, lat2, 'linear', NaN);
    
    if isnan(depth1) || isnan(depth2)
        warning('Impossible d''interpoler la profondeur aux points de contrôle. Méthode de repli utilisée.');
        [lat_interp, lon_interp] = fallback_interpolation(params.fallback, lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
        depth_path = [];
        return;
    end
    
    params.target_depth = (depth1 + depth2) / 2;
    %params.target_depth = -1888;
    
end

fprintf('Interpolation bathymétrique: méthode %s, profondeur cible %.0f m\n', params.method, params.target_depth);

%% Calcul du chemin selon la méthode choisie
switch params.method
    case 'isobath'
        [waypoints_lat, waypoints_lon] = find_isobath_path(lat1, lon1, lat2, lon2, lat_grid, lon_grid, depth_grid, params);
        
    otherwise
        error('Méthode inconnue: %s', params.method);
end

%% Vérification que le chemin a été trouvé
if isempty(waypoints_lat) || length(waypoints_lat) < 2
    warning('Impossible de trouver un chemin bathymétrique. Méthode de repli utilisée.');
    [lat_interp, lon_interp] = fallback_interpolation(params.fallback, lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
    depth_path = [];
    return;
end

%% Interpolation le long du chemin trouvé
[lat_interp, lon_interp] = interpolate_along_waypoints(waypoints_lat, waypoints_lon, juld_query, juld_start, juld_end);

%% Calcul des profondeurs le long du chemin (si demandé)
if nargout > 2
    depth_path = zeros(size(lat_interp));
    for i = 1:length(lat_interp)
        depth_path(i) = interp2(lon_grid, lat_grid, depth_grid, lon_interp(i), lat_interp(i), 'linear', NaN);
    end
end

end

%% ========================================================================
%% MÉTHODE 1: CHEMIN SUIVANT UNE ISOBATHE
%% ========================================================================
function [waypoints_lat, waypoints_lon] = find_isobath_path(lat1, lon1, lat2, lon2, lat_grid, lon_grid, depth_grid, params)
%FIND_ISOBATH_PATH Trouve un chemin suivant une isobathe

fprintf('  Recherche d''un chemin suivant l''isobathe %.0f m...\n', params.target_depth);

% Créer une grille plus fine pour la recherche de contours
[lon_fine, lat_fine] = meshgrid(linspace(min(lon_grid(:)), max(lon_grid(:)), size(lon_grid, 2)*2), ...
    linspace(min(lat_grid(:)), max(lat_grid(:)), size(lat_grid, 1)*2));
depth_fine = interp2(lon_grid, lat_grid, depth_grid, lon_fine, lat_fine, 'spline');

% Extraction des contours à la profondeur cible
target_depths = [params.target_depth - params.depth_tolerance, ...
    params.target_depth, ...
    params.target_depth + params.depth_tolerance];

best_path_lat = [];
best_path_lon = [];
min_distance = inf;

for target_dep = target_depths
    try
        % Extraction du contour
        figure(2)
        C = contour(lon_fine, lat_fine, depth_fine, [target_dep target_dep]);
        
        if size(C, 2) < 3
            continue;
        end
        
        % Analyser les contours trouvés
        [contour_segments] = parse_contour_matrix(C);
        
        % Trouver le segment le plus proche des points de départ/arrivée
        for i = 1:length(contour_segments)
            segment_lat = contour_segments{i}(:, 2);
            segment_lon = contour_segments{i}(:, 1);
            
            if length(segment_lat) < 3
                continue;
            end
            
            % Calculer la distance aux points de départ et d'arrivée
            dist_start = min(sqrt((segment_lat - lat1).^2 + (segment_lon - lon1).^2));
            dist_end = min(sqrt((segment_lat - lat2).^2 + (segment_lon - lon2).^2));
            total_dist = dist_start + dist_end;
            
            if total_dist < min_distance
                min_distance = total_dist;
                best_path_lat = segment_lat;
                best_path_lon = segment_lon;
            end
        end
        
    catch ME
        fprintf('    Erreur lors de l''extraction du contour %.0f m: %s\n', target_dep, ME.message);
        continue;
    end
end

% Si un contour a été trouvé, créer les points de passage
if ~isempty(best_path_lat)
    % Connecter le point de départ au contour le plus proche
    [~, idx_start] = min(sqrt((best_path_lat - lat1).^2 + (best_path_lon - lon1).^2));
    [~, idx_end] = min(sqrt((best_path_lat - lat2).^2 + (best_path_lon - lon2).^2));
    
    if idx_start >= idx_end
        best_path_lat=flipud(best_path_lat);
        best_path_lon=flipud(best_path_lon);
    end
    
    [~, idx_start] = min(sqrt((best_path_lat - lat1).^2 + (best_path_lon - lon1).^2));
    [~, idx_end] = min(sqrt((best_path_lat - lat2).^2 + (best_path_lon - lon2).^2));
    
    % Créer le chemin en suivant le contour
    %     if idx_start <= idx_end
    contour_lat = best_path_lat(idx_start:idx_end);
    contour_lon = best_path_lon(idx_start:idx_end);
    %     else
    %         contour_lat = best_path_lat([idx_start:end, 1:idx_end]);
    %         contour_lon = best_path_lon([idx_start:end, 1:idx_end]);
    %     end
    
    % Assembler le chemin complet
    waypoints_lat = [lat1; contour_lat; lat2];
    waypoints_lon = [lon1; contour_lon; lon2];
    
    % Ré-échantillonner pour avoir le nombre de waypoints souhaité
    if length(waypoints_lat) > params.n_waypoints
        indices = round(linspace(1, length(waypoints_lat), params.n_waypoints));
        waypoints_lat = waypoints_lat(indices);
        waypoints_lon = waypoints_lon(indices);
    end
    
    fprintf('    Chemin trouvé avec %d points de passage\n', length(waypoints_lat));
else
    fprintf('    Aucun chemin isobathique trouvé\n');
    waypoints_lat = [];
    waypoints_lon = [];
end

end


%% ========================================================================
%% FONCTIONS AUXILIAIRES
%% ========================================================================
function [lat_grid, lon_grid, depth_grid] = get_bathymetry_grid(lat1, lon1, lat2, lon2)
%GET_BATHYMETRY_GRID Récupère la grille bathymétrique avec m_tbase
marge=15;
% Définir les limites de la zone avec une marge
lat_min = min(lat1, lat2) - marge;
lat_max = max(lat1, lat2) + marge;
lon_min = min(lon1, lon2) - marge;
lon_max = max(lon1, lon2) + marge;

% Gestion de l'antimeridien
if abs(lon2 - lon1) > 180
    if lon1 > 0
        lon_min = min(lon1, lon2 + 360) - marge;
        lon_max = max(lon1, lon2 + 360) + marge;
    else
        lon_min = min(lon1 + 360, lon2) - marge;
        lon_max = max(lon1 + 360, lon2) + marge;
    end
end

% Contraindre aux limites valides
lat_min = max(-90, lat_min);
lat_max = min(90, lat_max);
lon_min = max(-180, lon_min);
lon_max = min(180, lon_max);

try
    % Créer la grille
    %     resolution = 0.1; % Résolution en degrés
    %     lon_vec = lon_min:resolution:lon_max;
    %     lat_vec = lat_min:resolution:lat_max;
    %     [lon_grid, lat_grid] = meshgrid(lon_vec, lat_vec);
    %
    % Récupérer la bathymétrie avec m_tbase
    fprintf('  Récupération de la bathymétrie (%.1f°x%.1f°)...\n', ...
        lon_max-lon_min, lat_max-lat_min);
    m_proj( 'mercator','longitude',[lon_min lon_max],'latitude',[lat_min lat_max]);
    %depth_grid = m_tbase(lon_grid, lat_grid);
    [depth_grid,lon_grid,lat_grid]=m_tbase([lon_min lon_max lat_min lat_max]);
    %[lon_grid, lat_grid] = meshgrid(lon_vec, lat_vec);
    % Vérifier que nous avons des données valides
    if all(isnan(depth_grid(:)))
        warning('Toutes les données bathymétriques sont NaN');
        depth_grid = [];
        return;
    end
    
    % Interpoler les NaN si nécessaire
    nan_mask = isnan(depth_grid);
    if any(nan_mask(:))
        fprintf('    Interpolation de %.1f%% de valeurs NaN\n', 100*sum(nan_mask(:))/numel(depth_grid));
        depth_grid = inpaint_nans(depth_grid, 4); % Méthode simple d'interpolation
    end
    
catch ME
    fprintf('  Erreur lors de la récupération de la bathymétrie: %s\n', ME.message);
    lat_grid = [];
    lon_grid = [];
    depth_grid = [];
end

end

function [contour_segments] = parse_contour_matrix(C)
%PARSE_CONTOUR_MATRIX Analyse la matrice de contours de MATLAB

contour_segments = {};
i = 1;
segment_idx = 1;

while i < size(C, 2)
    level = C(1, i);
    n_points = C(2, i);
    
    if i + n_points > size(C, 2)
        break;
    end
    
    % Extraire les points du segment
    segment_points = C(:, i+1:i+n_points)';
    contour_segments{segment_idx} = segment_points;
    
    i = i + n_points + 1;
    segment_idx = segment_idx + 1;
end

end

function [lat_interp, lon_interp] = interpolate_along_waypoints(waypoints_lat, waypoints_lon, juld_query, juld_start, juld_end)
%INTERPOLATE_ALONG_WAYPOINTS Interpole le long des points de passage

% Calculer les distances cumulées le long du chemin
distances = [0];
for i = 2:length(waypoints_lat)
    dist = sqrt((waypoints_lat(i) - waypoints_lat(i-1))^2 + (waypoints_lon(i) - waypoints_lon(i-1))^2);
    distances(i) = distances(i-1) + dist;
end

% Normaliser les distances
distances_norm = distances / max(distances);

% Facteur d'interpolation temporelle
t = (juld_query - juld_start) / (juld_end - juld_start);
t = max(0, min(1, t)); % Clamp entre 0 et 1

% Interpolation le long du chemin
lat_interp = interp1(distances_norm, waypoints_lat, t, 'pchip', 'extrap');
lon_interp = interp1(distances_norm, waypoints_lon, t, 'pchip', 'extrap');

end

function [lat_interp, lon_interp] = fallback_interpolation(method, lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end)
%FALLBACK_INTERPOLATION Méthodes de repli

switch method
    case 'linear'
        [lat_interp, lon_interp] = interpolate_linear_fallback(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
    case 'sphere'
        [lat_interp, lon_interp] = interpolate_sphere_fallback(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
    otherwise
        error('Méthode de repli inconnue: %s', method);
end

end

function [lat_interp, lon_interp] = interpolate_linear_fallback(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end)
%INTERPOLATE_LINEAR_FALLBACK Interpolation linéaire simple

f = (juld_query - juld_start) / (juld_end - juld_start);

% Gestion de l'antimeridien
dlon = lon2 - lon1;
if abs(dlon) > 180
    if dlon > 0
        lon1 = lon1 + 360;
    else
        lon2 = lon2 + 360;
    end
end

lat_interp = lat1 + f .* (lat2 - lat1);
lon_interp = lon1 + f .* (lon2 - lon1);

% Normalisation des longitudes
lon_interp = mod(lon_interp + 180, 360) - 180;

end

function [lat_interp, lon_interp] = interpolate_sphere_fallback(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end)
%INTERPOLATE_SPHERE_FALLBACK Interpolation sphérique simple

% Conversion en radians
lat1_rad = deg2rad(lat1);
lon1_rad = deg2rad(lon1);
lat2_rad = deg2rad(lat2);
lon2_rad = deg2rad(lon2);

% Vecteurs 3D
p1 = [cos(lat1_rad)*cos(lon1_rad), cos(lat1_rad)*sin(lon1_rad), sin(lat1_rad)];
p2 = [cos(lat2_rad)*cos(lon2_rad), cos(lat2_rad)*sin(lon2_rad), sin(lat2_rad)];

% Angle entre les points
omega = acos(max(-1, min(1, dot(p1, p2))));

% Facteur d'interpolation
f = (juld_query - juld_start) / (juld_end - juld_start);

% Interpolation SLERP
if omega < 1e-10
    lat_interp = repmat(lat1, size(f));
    lon_interp = repmat(lon1, size(f));
else
    sin_omega = sin(omega);
    A = sin((1 - f) * omega) / sin_omega;
    B = sin(f * omega) / sin_omega;
    
    p_interp = A(:) * p1 + B(:) * p2;
    
    % Conversion en lat/lon
    x = p_interp(:, 1);
    y = p_interp(:, 2);
    z = max(-1, min(1, p_interp(:, 3)));
    
    lat_interp = rad2deg(asin(z))';
    lon_interp = rad2deg(atan2(y, x))';
end

% Normalisation des longitudes
lon_interp = mod(lon_interp + 180, 360) - 180;

end

%% ========================================================================
%% FONCTION D'INTERPOLATION NaN SIMPLE
%% ========================================================================
function B = inpaint_nans(A, method)
%INPAINT_NANS Interpolation simple des NaN (version basique)

if nargin < 2
    method = 4; % Méthode par défaut
end

B = A;
[m, n] = size(A);

% Masque des NaN
nanmask = isnan(A);

if ~any(nanmask(:))
    return; % Aucun NaN à interpoler
end

% Méthode simple: remplacer par la moyenne des voisins non-NaN
for i = 1:m
    for j = 1:n
        if nanmask(i, j)
            % Définir le voisinage
            i_min = max(1, i-1);
            i_max = min(m, i+1);
            j_min = max(1, j-1);
            j_max = min(n, j+1);
            
            % Extraire les valeurs du voisinage
            neighbors = A(i_min:i_max, j_min:j_max);
            valid_neighbors = neighbors(~isnan(neighbors));
            
            if ~isempty(valid_neighbors)
                B(i, j) = mean(valid_neighbors);
            end
        end
    end
end

% Si il reste des NaN, utiliser la valeur médiane globale
remaining_nans = isnan(B);
if any(remaining_nans(:))
    global_median = nanmedian(A(:));
    if ~isnan(global_median)
        B(remaining_nans) = global_median;
    else
        B(remaining_nans) = -3000; % Profondeur par défaut
    end
end

end

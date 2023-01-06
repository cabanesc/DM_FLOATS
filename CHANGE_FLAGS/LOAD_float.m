% -========================================================
%   USAGE : LOAD_float(floatname,dacname,varargin)
%   PURPOSE : plot Argo profiles and modify quality flags in netcdf files using a GUI
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char) e.g.  'coriolis'
%
%   OPTIONNAL INPUT :
%    'ERASE' (logical)  ERASE=1 if float netcdf files (stored in  a local directory)  are replaced by those available on GDAC, ERASE = 0 (default) if not.
%    'ASK'   (logical) ASK=1 (default) if a question dialog bow is open to confirm that you want to reload the file
%                          0  NO DIALOG BOX
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created (2016) ccabanes
%
%   CALLED SUBROUTINES:
% -------------------------------------
% GIT BRANCH: PourOXY
% ========================================================
function LOAD_float(floatname,dacname,varargin)


%init_path; % Commented by Reynaud
n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);
ERASE=0;
ASK=1;

if isfield(s,'ERASE')==1;ERASE=s.ERASE;end;
if isfield(s,'ASK')==1;ASK=s.ASK;end;


C=load_configuration('config.txt');

DIR_FTP_CORIOLIS=C.DIR_FTP_CORIOLIS;
CONFIG.FILE_TOPO=C.FILE_TOPO;
CONFIG.DIR_FTP=C.DIR_FTP;


% FLOTEUR ANALYSE
CONFIG.floatname=floatname;
CONFIG.dacname=dacname;
%  n=length(varargin);
%  
%  if n/2~=floor(n/2)
%      error('check the imput arguments')
%  end
%  
%  f=varargin(1:2:end);
%  c=varargin(2:2:end);
%  s = cell2struct(c,f,2);
%  ERASE=0;
%  VPN=1;
%  FLAG=1;
%  if isfield(s,'ERASE')==1;ERASE=s.ERASE;end;
%  if isfield(s,'VPN')==1;VPN=s.VPN;end;
%  if isfield(s,'FLAG')==1;FLAG=s.FLAG;end;
%  
%  CONFIG.VPN=VPN;
%  CONFIG.FLAG=FLAG;

if ERASE==1 & exist([CONFIG.DIR_FTP '/' CONFIG.dacname '/' CONFIG.floatname])
    if ASK==1
    isw = questdlg('Do you want to reload the float netcdf files ?','Loading float files...','YES', 'NO','NO');
    else
    isw='YES';
    end
    switch isw
        case {'YES'} 
    destinf=[CONFIG.DIR_FTP '/' CONFIG.dacname '/' CONFIG.floatname '/'];
    
    status0=rmdir(destinf,'s');
    end
end

if ~exist([CONFIG.DIR_FTP '/' CONFIG.dacname '/' CONFIG.floatname])
    sourcef=[DIR_FTP_CORIOLIS '/' CONFIG.dacname '/' CONFIG.floatname '/*.nc'];
    destinf=[CONFIG.DIR_FTP '/' CONFIG.dacname '/' CONFIG.floatname '/'];
    status1=copyfile(sourcef,destinf);
    disp(['copy status (1 is ok): ' num2str(status1) ])
end
if ~exist([CONFIG.DIR_FTP '/' CONFIG.dacname '/' CONFIG.floatname '/profiles/'])
    sourcef=[DIR_FTP_CORIOLIS '/' CONFIG.dacname '/' CONFIG.floatname '/profiles/R*.nc'];
    destinf=[CONFIG.DIR_FTP '/' CONFIG.dacname '/' CONFIG.floatname '/profiles/'];
    status2=copyfile(sourcef,destinf);
    sourcef=[DIR_FTP_CORIOLIS '/' CONFIG.dacname '/' CONFIG.floatname '/profiles/D*.nc'];
    destinf=[CONFIG.DIR_FTP '/' CONFIG.dacname '/' CONFIG.floatname '/profiles/'];
    status3=copyfile(sourcef,destinf);
    disp(['copy status Core R & D Files (1 is ok): '  num2str(status2) ', ' num2str(status3)])
	
    sourcef=[DIR_FTP_CORIOLIS '/' CONFIG.dacname '/' CONFIG.floatname '/profiles/BR*.nc'];
    destinf=[CONFIG.DIR_FTP '/' CONFIG.dacname '/' CONFIG.floatname '/profiles/'];
    status2=copyfile(sourcef,destinf);
    sourcef=[DIR_FTP_CORIOLIS '/' CONFIG.dacname '/' CONFIG.floatname '/profiles/BD*.nc'];
    destinf=[CONFIG.DIR_FTP '/' CONFIG.dacname '/' CONFIG.floatname '/profiles/'];
    status3=copyfile(sourcef,destinf);
    disp(['copy status BGC R & D Files(1 is ok): '  num2str(status2) ', ' num2str(status3)])
end

eval('cd ..');



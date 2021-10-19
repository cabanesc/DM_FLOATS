% FUNCTION [DATE_CORIO] = DATE_CONV(DATE_GREG, [FORM])
%
% Convert a gregorian date or array of dates into string according to the 
% format form 
% 
% Input parameters 
%	DATE_GREG : gregorian date 
%       date_greg(:,1) = year 4 digits
%       date_greg(:,2) = month
%       date_greg(:,3) = day
%       date_greg(:,4) = hour
%       date_greg(:,5) = minute
%       date_greg(:,6) = seconds
%
%	FORM : requested output format 
%       default      ==> 'DD/MM/YYYY HH24:MI:SS'
%  		'c[oriolis]'    ==> 'DD/MM/YYYY HH24:MI:SS'
%		'm[atlab]'      ==> 'DD-MMM-YYYY HH24:MM:SS' 
%		'a[rgo]'        ==> 'YYYYDDMMHH24MISS' 
%
%  Output parameters :
%  		date_greg_corio : date string in the requested format
%
%
% T.Terre    07/02/2001 adapted from date_conv.m used for CORIOLIS
%
% C. Lagadec 28/07/02   ajout du traitement de la date convention ARGO
%
% -----------------------------------------------------------------

function [date_corio] = date_conv(date_greg, varargin);

%
% Check input parameters - Force to default eventually
%
date_corio = [];
year = [];
day  = [];
hou  = [];
minu = [];
secd = [];

if (nargin == 1)
	form = 'C';
end

if (nargin == 2)
	if (ischar(varargin{1}))
		form = upper(varargin{1});
	else
		disp('Need a string as 2nd parameter : C[ORIOLIS] or M[ATLAB] or A[RGO]');
		return;
	end
end


%  tronque les secondes:
xsec  = round(date_greg(:,6));

year  = num2str(date_greg(:,1), '%4.4i');
day   = num2str(date_greg(:,3), '%2.2i');
hour  = num2str(date_greg(:,4), '%2.2i');
minu  = num2str(date_greg(:,5), '%2.2i');
secd  = num2str(xsec, '%2.2i');
 
%
% Select the correct format
%
[nn,m] = size(date_greg);


%
% Form return string
% 

dl_hour = ones(nn,1)*':';
dl_blan = ones(nn,1)*' ';


switch form(1)
	case 'C'
		dl_year = ones(nn,1)*'/';
		month = num2str(date_greg(:,2), '%2.2i');
                date_corio  = ...
 	   [day dl_year month dl_year year dl_blan hour dl_hour minu dl_hour secd];

	case 'M'
		dl_year = ones(nn,1)*'-';
		            month=datestr(datenum(date_greg(:,1),date_greg(:,2),date_greg(:,3)),3);
                date_corio  = ...
 	   [day dl_year month dl_year year dl_blan hour dl_hour minu dl_hour secd];

	case 'A'
		month = num2str(date_greg(:,2), '%2.2i');
                date_corio  =  [year month day hour minu secd];

	otherwise
		disp('       Unknown conversion format');
		return
end

	

% FUNCTION NUMDATE = DATESTR2NUM(DATE_TIME, [FORM])
%
% Convert a string date or array of string dates into gregorian dates or 
% a julian days starting and ending at midnight.
%
% Input parameters : 
% 		DATE_TIME : date in any of the formats below 
%					'jj/mm/yyyy hh:mm:ss'
%                   'jj-mmm-yyyy hh:mm:ss'
%
%		FORM      : Output format for the numerical value
%					'G[REGORIAN]' ==> numdate = [yyyy mo da hr mi sec]
%                   'J[ULIAN]     ==> numdate = Julian day start and end at 
%                                               midnight.
%
% Output parameter :
%		NUMDATE   : Converted value accordind to desired format 
% 
%    In this convention, Julian day 2440000 began at 0000 hours, May 23, 1968.
%

% T.Terre 07/02/2001 

function numdate = datestr2num(date_time, varargin)

%
% Check input parameters - Force to default eventually
%
numdate = [];
if (nargin == 1)
	form = 'GREG';
end

if (nargin == 2)
	if (ischar(varargin{1}))
		form = upper(varargin{1});
	else
		disp('Need a string as 2nd parameter : G[REGORIAN] or J[ULIAN]');
		return;
	end
end

%
% Determine length of array
%
nn = size(date_time, 1);

%
% Process according to input format 
% i.e detect if / or - is present twice
%	   'jj/mm/yyyy hh:mm:ss'
%   or 'JJ-MMM-YYYY hh:mm:ss'
%
dt = zeros(nn, 8);
dtime = zeros(nn, 6);
for i=1:nn
	if (length(findstr(date_time(i, :), '/')) == 2)
		a = sscanf(date_time(i,:),'%d/%d/%d %d:%d:%d')';
		dt(i, 1:length(a)) = a;
		dtime(i,:) = [dt(i, 3) dt(i, 2) dt(i, 1) dt(i, 4), dt(i, 5) dt(i, 6)];
	else if (length(findstr(date_time(i, :), '-')) == 2)
			a = sscanf(date_time(i, :), '%d-%3c-%d %d:%d:%d')';
			dt(i, 1:length(a)) = a;  
			month = upper(char(dt(i, 2:4)));
			switch month
				case 'JAN'
					mm = 1;
				case 'FEB'
					mm = 2;
				case 'MAR'
					mm = 3;
				case 'APR'
					mm = 4;
				case 'MAY'
					mm = 5;
				case 'JUN'
					mm = 6;
				case 'JUL'
					mm = 7;
				case 'AUG'
					mm = 8;
				case 'SEP'
					mm = 9;
				case 'OCT'
					mm = 10;
				case 'NOV'
					mm = 11;
				case 'DEC'
					mm = 12;
			end
			dtime(i,:) = [dt(i, 5) mm dt(i, 1) dt(i, 6) dt(i, 7) dt(i, 8)];
		else
			disp('Input string not correct : jj/mm/yyyy hh:mm:ss or JJ-MMM-YYYY hh:mm:ss')
			return
		end
	end
end

%
% Form the return value 
%
switch form(1)
	case 'G'
		numdate = greg_0h(jul_0h(dtime));

	case 'J'
		numdate = jul_0h(dtime);

	otherwise
		disp('       Unknown conversion format');
end

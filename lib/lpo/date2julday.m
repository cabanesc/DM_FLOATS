% FUNCTION JUL_DAY = DATE2JULDAY(DATE_TIME, [FORM])
%
% Convert a string date or array of string dates given in format indicated by 
% form into a julian day starting and ending at midnight.
%
% Input parameter : 
% 		date_time : date in format described by the parameter form
%		form      : none ==> 'jj/mm/yyyy hh:mm:ss'
%					0    ==> 'jj/mm/yyyy hh:mm:ss'
%                   1    ==> 'JJ-MMM-YYYY hh:mm:ss'
%
% Output parameter :
%		jul_day : julian day referenced to hour 0.
% 
%    In this convention, Julian day 2440000 began at 0000 hours, May 23, 1968.
%

% T.Terre 07/02/2001 

function jul_day = date2julday(date_time, varargin)

if (nargin == 1)
	form = 0;
end

if (nargin == 2)
	form = varargin{1};
end

nn = size(date_time, 1);

switch form
	case 0
		dt = zeros(nn, 6);
		for i=1:nn
			dt(i, :) = sscanf(date_time(i,:),'%d/%d/%d %d:%d:%d')';
		end
		dtime = [dt(:, 3) dt(:, 2) dt(:, 1) dt(:, 4), dt(:, 5) dt(:, 6)];

	case 1
		dt = zeros(nn, 8);
		mm = zeros(nn,1);
		for i=1:nn
			dt(i, :)  = sscanf(date_time(i, :), '%d-%3c-%d %d:%d:%d')';

			month = upper(char(dt(i, 2:4)));
			switch month
				case 'JAN'
				mm(i) = 1;

				case 'FEB'
				mm(i) = 2;

				case 'MAR'
				mm(i) = 3;

				case 'APR'
				mm(i) = 4;

				case 'MAY'
				mm(i) = 5;

				case 'JUN'
				mm(i) = 6;

				case 'JUL'
				mm(i) = 7;

				case 'AUG'
				mm(i) = 8;

				case 'SEP'
				mm(i) = 9;

				case 'OCT'
				mm(i) = 10;

				case 'NOV'
				mm(i) = 11;

				case 'DEC'
				mm(i) = 12;

			end
			dtime = [dt(:, 5) mm dt(:, 1) dt(:, 6) dt(:, 7) dt(:, 8)];
		end
	end
jul_day = jul_0h(dtime);

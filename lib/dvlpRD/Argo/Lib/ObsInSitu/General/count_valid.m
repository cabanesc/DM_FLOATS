function Test=count_valid(Co,onechamp)
% -========================================================
%   USAGE : [Test]=check_isfillval_prof(Co)
%   PURPOSE : short description of the function here
% -----------------------------------
%   INPUT :
%     IN1   (class)  -comments-
%             additional description
%     IN2   (class)  -comments-
%
%   OPTIONNAL INPUT :
%    OPTION1  (class)  -comments-
% -----------------------------------
%   OUTPUT :
%     OUT1   (class)  -comments-
%             additional description
%     OUT2   (class)  -comments-
%             additional description
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx 
%   CALLED SUBROUTINES: none
% ========================================================


if isempty(strfind(Co.obj,'ObsInSitu'))
   error('count_valid not define for this type of structure')
else
    Test=[];
    
    
	if isfield(Co.(onechamp),'data')
	    if isempty(Co.(onechamp).data)==0
		Test= logical(zeros(size(Co.(onechamp).data,1),1));
		if isempty(Co.(onechamp).FillValue_)==0
		    Test= sum(Co.(onechamp).data~=Co.(onechamp).FillValue_,2);
		end
	    end
	end
end
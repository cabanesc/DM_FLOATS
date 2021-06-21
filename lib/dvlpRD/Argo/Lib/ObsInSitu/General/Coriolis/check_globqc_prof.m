function [Testi,GlobQC]=check_globqc_prof(Co,onechamp,fillval)
% -========================================================
%   USAGE : [Test]=check_isfillval_prof(Co)
%           [Testi,GlobQC]=check_globqc_prof(Co,onechamp,fillval)
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

if isempty(strfind(Co.obj,'ObsInSitu/Coriolis'))
    %error('check_globqc_prof not define for this type of structure')
else
    Testi=[];
    GlobQC=[];
    if nargin==2
        fillval='FillValue_';
    end
    if isfield(Co,onechamp)
        if isfield(Co.(onechamp),'data')
            if isempty(Co.(onechamp).data)==0
                qc = Co.(onechamp).data;
                if isempty(Co.(onechamp).(fillval))==0

                    isfill=(qc==Co.(onechamp).(fillval)|qc==9);
                    qc(qc==1|qc==2|qc==5|qc==8)=1;
                    qc(qc==3|qc==4|qc==0)=2;
                    A=sum(qc==1,2);
                    B=sum(~isfill,2);
                    nf=B~=0;
                    Testi=999*ones(size(qc,1),1);
                    GlobQC=repmat(' ',[size(qc,1),1]);
                    
                    Testi(nf)=A(nf)./B(nf);
                    GlobQC(Testi==0)='F';
                    GlobQC(Testi>0)='E';
                    GlobQC(Testi>=0.25)='D';
                    GlobQC(Testi>=0.5)='C';
                    GlobQC(Testi>=0.75)='B';
                    GlobQC(Testi==1)='A';
                    GlobQC(Testi==999)=' ';
                    %  	        Testi(T==0)=6;
                    %  	        Testi(T>0)=5;
                    %  	        Testi(T>0.25)=4;
                    %  	        Testi(T>0.5)=3;
                    %  	        Testi(T>0.75)=2;
                    %  	        Testi(T==1)=1;
                    %  	        Testi(T==9)=0;
                end
            end
        end
    end
end



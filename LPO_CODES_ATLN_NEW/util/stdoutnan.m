function y = std(x,flag,dim)
%STD    Standard deviation.
%   For vectors, STD(X) returns the standard deviation. For matrices,
%   STD(X) is a row vector containing the standard deviation of each
%   column.  For N-D arrays, STD(X) is the standard deviation of the
%   elements along the first non-singleton dimension of X.
%
%   STD(X) normalizes by (N-1) where N is the sequence length.  This
%   makes STD(X).^2 the best unbiased estimate of the variance if X
%   is a sample from a normal distribution.
%
%   STD(X,1) normalizes by N and produces the second moment of the
%   sample about its mean.  STD(X,0) is the same as STD(X).
%
%   STD(X,FLAG,DIM) takes the standard deviation along the dimension
%   DIM of X.  When FLAG=0 STD normalizes by (N-1), otherwise STD
%   normalizes by N.
%
%   Example: If X = [4 -2 1
%                    9  5 7]
%     then std(X,0,1) is [ 3.5355 4.9497 4.2426] and std(X,0,2) is [3.0
%                                                                   2.0]
%   See also COV, MEAN, MEDIAN, CORRCOEF.

%   J.N. Little 4-21-85
%   Revised 5-9-88 JNL, 3-11-94 BAJ, 5-26-95 dlc, 5-29-96 CMT.
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.17 $  $Date: 1997/11/21 23:24:08 $

if nargin<2, flag = 0; end
if nargin<3, 
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end

% Avoid divide by zero.
if size(x,dim)==1, y = zeros(size(x)); return, end

tile = ones(1,max(ndims(x),dim));
tile(dim) = size(x,dim);

%xc = x - repmat(sum(x,dim)/size(x,dim),tile);  % Remove mean
xc = x - repmat(meanoutnan(x,dim),tile);  % Remove mean
xcnan=(isnan(xc));
coef=sum(~xcnan,dim);
coef(coef==0)=NaN;
xc(xcnan)=0;
if flag,
  y = sqrt(sum(conj(xc).*xc,dim)./coef);
else
  coef(coef==1)=2; %then y=0
  y = sqrt(sum(conj(xc).*xc,dim)./(coef-1));
end


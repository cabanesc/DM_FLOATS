function [zf,thf]=idfilt(z,n,Wn,hs)
%IDFILT Filters data.
%   ZF = IDFILT(Z,N,Wn)
%
%   Z  : The output-input data matrix Z = [Y U]. Multi-variable
%        data are allowed
%   ZF : The filtered data. The columns of ZF correspond to those of Z.
%   N  : The order of the filter to be used.
%        N = a positive integer: gives an acausal, zero-phase filter of
%            order 2N, determined from a Butterworth filter.
%        N = a negative integer: gives a causal Butterworth filter of
%            order abs(N).
%   Wn : The cut-off frequency(ies) in fractions of the Nyquist frequency.
%        If Wn is a scalar a Low-Pass filter will be used.
%        If Wn = [Wl Wh] a Band-Pass filter with pass-band  Wn will be used
%
%   With ZF = IDFILT(Z,N,Wn,'high') or ZF = IDFILT(Z,N,[Wl Wh],'stop')
%   High-Pass and Band-Stop filters will be used instead.
%
%   With [ZF,THFILT] = IDFILT(Z,N,Wn), also the filter will be returned
%   in the THETA-format as THFILT. (See also THETA).
%   See also DTREND and IDRESAMP.

if nargin<3
   disp('Usage: zf = IDFILT(DATA,ORDER,BAND);')
   disp('       zf = IDFILT(DATA,ORDER,BAND,''HIGH'');')
   disp('       with ORDER = an integer (positive for acausal filtering).')
   return
end
if nargin<4,hs=[];end
if length(Wn)>1,
   if Wn(1)>=Wn(2)
     error('The band is ill-defined. Should be Wn = [a b] with a < b.')
   end
   if abs(Wn(1))<eps
      Wn=Wn(2);
   elseif abs(Wn(2)-1)<eps,
      Wn=Wn(1);
      if isempty(hs),hs='high';else hs=[];end
   end
end

[ndat,nyu]=size(z);
if ndat<nyu
   error('Data should be organized columnwise.')
end
if ndat<3*n
 error('The data record length must be at least three times the filter order')
end
if n<0,n=abs(n);causal=1;else causal=0;end

btype = 1;
if ~isempty(hs),btype=3;end

if length(Wn) == 2
        btype = btype + 1;
end

% step 1: get analog, pre-warped frequencies
        fs = 2;
        u = 2*fs*tan(pi*Wn/fs);

% step 2: convert to low-pass prototype estimate
if btype == 1   % lowpass
        Wn = u;
elseif btype == 2       % bandpass
        Bw = u(2) - u(1);
        Wn = sqrt(u(1)*u(2));   % center frequency
elseif btype == 3       % highpass
        Wn = u;
elseif btype == 4       % bandstop
        Bw = u(2) - u(1);
        Wn = sqrt(u(1)*u(2));   % center frequency
end

% step 3: Get N-th order Butterworth analog lowpass prototype
p = exp(sqrt(-1)*(pi*(1:2:2*n-1)/(2*n) + pi/2)).';
k = real(prod(-p));

% Transform to state-space
% Strip infinities and throw away.
p = p(finite(p));

% Group into complex pairs
np = length(p);
nz= 0;%= length(z);
p = cplxpair(p,1e6*np*norm(p)*eps + eps);

% Initialize state-space matrices for running series
a=[]; b=[]; c=zeros(1,0); d=1;

% If odd number of poles only, convert the pole at the
% end into state-space.
%  H(s) = 1/(s-p1) = 1/(s + den(2))
if rem(np,2)
        a = p(np);
        b = 1;
        c = 1;
        d = 0;
        np = np - 1;
end

% Now we have an even number of poles and zeros, although not
% necessarily the same number - there may be more poles%.
%   H(s) = (s^2+num(2)s+num(3))/(s^2+den(2)s+den(3))
% Loop thru rest of pairs, connecting in series to build the model.
i = 1;

% Take care of any left over unmatched pole pairs.
%   H(s) = 1/(s^2+den(2)s+den(3))
while i < np
        den = real(poly(p(i:i+1)));
        wns = sqrt(prod(abs(p(i:i+1))));
        if wns == 0, wns = 1; end
        t = diag([1 1/wns]);    % Balancing transformation
        a1 = t\[-den(2) -den(3); 1 0]*t;
        b1 = t\[1; 0];
        c1 = [0 1]*t;
        d1 = 0;
%       [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1);
% Next lines perform series connection
        [ma1,na1] = size(a);
        [ma2,na2] = size(a1);
        a = [a zeros(ma1,na2); b1*c a1];
        b = [b; b1*d];
        c = [d1*c c1];
        d = d1*d;

        i = i + 2;
end
% Apply gain k:
c = c*k;
d = d*k;
% step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn
wo = Wn;
if btype == 1           % Lowpass
        at = wo*a;
        bt = wo*b;
        ct = c;
        dt = d;

elseif btype == 2       % Bandpass
        [ma,nb] = size(b);
        [mc,ma] = size(c);

        % Transform lowpass to bandpass
        q = wo/Bw;
        at = wo*[a/q eye(ma); -eye(ma) zeros(ma)];
        bt = wo*[b/q; zeros(ma,nb)];
        ct = [c zeros(mc,ma)];
        dt = d;

elseif btype == 3       % Highpass

        at =  wo*inv(a);
        bt = -wo*(a\b);
        ct = c/a;
        dt = d - c/a*b;

elseif btype == 4       % Bandstop
        [ma,nb] = size(b);
        [mc,ma] = size(c);

        % Transform lowpass to bandstop
        q = wo/Bw;
        at =  [wo/q*inv(a) wo*eye(ma); -wo*eye(ma) zeros(ma)];
        bt = -[wo/q*(a\b); zeros(ma,nb)];
        ct = [c/a zeros(mc,ma)];
        dt = d - c/a*b;

end
a=at;b=bt;c=ct;d=dt;
% step 5: Use Bilinear transformation to find discrete equivalent:

        t = 1/fs;
        r = sqrt(t);
        t1 = eye(size(a)) + a*t/2;
        t2 = eye(size(a)) - a*t/2;
        ad = t2\t1;
        bd = t/r*(t2\b);
        cd = r*c/t2;
        dd = c/t2*b*t/2 + d;
        a = ad;b = bd; c = cd; d = dd;
        den = poly(a);
        num = poly(a-b*c)+(d-1)*den;
   if nargout>1, thf=mktheta(den,num);end
   if causal
       for k=1:nyu
          x=ltitr(a,b,z(:,k));
          zf(:,k)=x*c'+d*z(:,k);
       end
   else
      for k = 1:nyu
             x=z(:,k);
             len = size(x,1);   % length of input
             b = num(:).';
             a = den(:).';
             nb = length(b);
             na = length(a);
             nfilt = max(nb,na);

             nfact = 3*(nfilt-1);  % length of edge transients

             if nb < nfilt, b(nfilt)=0; end   % zero-pad if necessary
             if na < nfilt, a(nfilt)=0; end
             zi = ( eye(nfilt-1) - [-a(2:nfilt).' ...
                  [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ...
                 ( b(2:nfilt).' - a(2:nfilt).'*b(1) );

             y = [2*x(1)-x((nfact+1):-1:2);x;...
                      2*x(len)-x((len-1):-1:len-nfact)];

            % filter, reverse data, filter again, and reverse data again
            y = filter(b,a,y,[zi*y(1)]);
            y = y(length(y):-1:1);
            y = filter(b,a,y,[zi*y(1)]);
            y = y(length(y):-1:1);

            % remove extrapolated pieces of y
            y([1:nfact len+nfact+(1:nfact)]) = [];
            zf(:,k)=y;
      end
   end

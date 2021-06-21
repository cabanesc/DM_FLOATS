% fonction qui calcule les fonctions d'hermite
% meridiens
% adaptation matlab du code fvalue2.f

function [fnh,fnu]=calpsi(ynorth,cphi,nyi);

nmodm=10;
nmp1=nmodm+1;
nmp2=nmodm+2;

kmprdg=111000.;
pi=4.*atan(1.);
beta=2.28671e-11;
xscale=sqrt(cphi/(beta));

ybot=-ynorth*kmprdg;
dy=2.*ynorth/(nyi-1);
hx=dy*kmprdg;
for j=1:nyi
  y(j)=ybot+(j-1)*hx;
end


% calculate hermite functions psi's.  Needs to be calculate
% for both vertical modes because of non-dimensional y.     
      
for ihmode=1:nmp2
    ihm=ihmode-1;

% calculate factorial associated with ihmode
    xf=1.;
    if(ihm~=0)
      for ih=1:ihm
        xf=xf*ih;
      end
    end

    for iy=1:nyi

% calculate nondimensional y
      ynd=y(iy)/xscale;

% calculate appropiate psi terms

      factor1=1/sqrt(sqrt(pi) * (2^ihm) * xf);
      factor1=factor1*exp(-(ynd*ynd)/2.);
      ih1=ihmode;
      if(ih1==1) 
        xh=1.;
      end
      if(ih1==2) 
        xh=2.*ynd;
      end
      if(ih1==3) 
        xh=(4.*(ynd^2))-2.;
      end
      if(ih1==4) 
        xh=(8.*(ynd^3))-(12.*ynd);
      end
      if(ih1==5) 
        xh=(16.*(ynd^4))-(48.*(ynd^2))+12.;
      end
      if(ih1==6) 
        xh=(32.*(ynd^5))-(160.*(ynd^3))+(120.*ynd);
      end
      if(ih1==7) 
        xh=(64. *(ynd^6))-(480.*(ynd^4))+(720.*(ynd^2))-120.;
      end
      if(ih1==8) 
        xh= 128.*(ynd^7)-1344.*(ynd^5)+3360.*(ynd^3)-1680.* ynd;
      end
      if(ih1==9) 
        xh=  256.*(ynd^8)-3584.*(ynd^6)+13440.*(ynd^4)-13440.*(ynd^2)+1680.;
      end
      if(ih1==10) 
        xh=  512.*(ynd^9)- 9216.*(ynd^7)+48384.*(ynd^5)-80640.*(ynd^3)+30240.* ynd;
      end
      if(ih1==11) 
        xh=  1024.*(ynd^10)-23040.*(ynd^8)+161280.*(ynd^6) -403200.*(ynd^4)+302400.*(ynd^2)-30240.;
      end
      if(ih1==12) 
        xh=   2048.*(ynd^11)- 56320.*(ynd^9)+506880.*(ynd^7) -1774080.*(ynd^5)+2217600.*(ynd^3)-665280.* ynd;
      end
      if(ih1==13) 
        xh=   4096.*(ynd^12)- 135168.*(ynd^10)+1520640.*(ynd^8)-7096320.*(ynd^6)+13305600.*(ynd^4)- 7983360.*(ynd^2) +665280.;
      end 
      psi(ihmode,iy)=factor1*xh;
    end
end

% start foring calculation of h for this particular x and t over
% all y grids

for iy=1:nyi

% for  kelvin wave terms

  fnh(1,iy)=(psi(1,iy)/sqrt(2.));
  fnu(1,iy)=(psi(1,iy)/sqrt(2.));
  fnv(1,iy)=0.0;

% now  rossby wave terms
 
  imo=1;
  for ih=1:nmp2-2
    imo=imo+1;
    arg=2.;
    rf=-0.5;

    factor2=((ih+1)^(rf))*psi(ih+2,iy);
    factor3= ( ih^(rf)) *psi(ih,iy);

% factor for v
    
    factor4=psi(ih+1,iy);
    
    fnu(imo,iy)=(factor2-factor3)/(2.*sqrt(arg));
    fnh(imo,iy)=(factor2+factor3)/(2.*sqrt(arg)); 
    fnv(imo,iy)=factor4;  
    
   end
end

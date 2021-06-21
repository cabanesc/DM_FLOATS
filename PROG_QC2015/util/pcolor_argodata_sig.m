function [hf] = pcolor_argodata_sig(cycnum,pres,sig0,param,paramname,vshading);
% Fonction qui trace une section Argo avec pcolor
% cycnum= tableau 1D des numeros de cycle
% pres = tableau 2D de pression
% param = tableau 2D du parametre a tracer en focntion de la pression
% paramname = nom du parametre:
%           = 'PSAL', 'TEMP', TPOT', 'DOXY' ou 'SIGO'
% vshading = definit le type de shading applique au pcolor
%          = 'faceted'
%          = 'interp'; un contour est aussi trace
%          = 'flat'


[nz,nprf]=size(pres);
% interpolation sur une grille pres reguliere.
regpress=repmat([1:1:max(max(pres))]',1,nprf);
for ik=1:nprf
    iin=~isnan(pres(:,ik));
    if sum(iin)>=2
    regparam(:,ik)=interp1(pres(iin,ik),param(iin,ik),regpress(:,ik));
    else
    regparam(:,ik)=NaN*regpress(:,ik);
    end
end
for ik=1:nprf
    iin=~isnan(pres(:,ik));
    if sum(iin)>=2
    regsig0(:,ik)=interp1(pres(iin,ik),sig0(iin,ik),regpress(:,ik));
    else
    regsig0(:,ik)=NaN*regpress(:,ik);
    end
end
param=regparam;
pres=regpress;
sig0=regsig0;
%param=regparam(1:10:end,:);
%pres=regpress(1:10:end,:);
[nz,nprf]=size(pres);
tabpres=NaN*ones(nz+1,nprf+1);
tabpres(1:nz,1:nprf)=pres;
tabpres(:,nprf+1)=tabpres(:,nprf);
tabpres(nz+1,:)=tabpres(nz,:);
tabparam=NaN*ones(nz+1,nprf+1);
tabparam(1:nz,1:nprf)=param;
tabparam(:,nprf+1)=tabparam(:,nprf);
tabparam(nz+1,:)=tabparam(nz,:);
tabsig=NaN*ones(nz+1,nprf+1);
tabsig(1:nz,1:nprf)=sig0;
tabsig(:,nprf+1)=tabsig(:,nprf);
tabsig(nz+1,:)=tabsig(nz,:);
    %keyboard


% defintion des intervalles
switch paramname
    case 'PSAL'
       coefin=10;
       coefout=1;
       pas=0.01;
       %pas=0.05;
    case {'TEMP' , 'TPOT'} 
        coefin=1;
        coefout=1;
        pas=0.2;
        %pas=0.5;
    case 'DOXY'
        coefin=1;
        coefout=10;
        pas=10;
    case 'OXYL'
        coefin=1;
        coefout=10;
        pas=10;
    case 'SIG0'
        coefin=10;
        coefout=1;
        pas=0.01; 
        %pas=0.1; 
    case 'PSAT'
        coefin=1;
        coefout=10;
        pas=5;
end
% paramname
% keyboard
cmin=coefout*min(floor(param(pres>1000)*coefin/coefout))/coefin;
cmax=coefout*max(ceil(param(pres>1000)*coefin/coefout))/coefin;
tabval2=[cmin-0.0001:pas:cmax+0.0001];

cmin=coefout*min(min(floor(param(pres>200)*coefin/coefout))/coefin);
cmax=coefout*max(max(ceil(param(pres>200)*coefin/coefout))/coefin);
%cmin=coefout*min(min(floor(param(pres>0)*coefin/coefout))/coefin);
%cmax=coefout*max(max(ceil(param(pres>0)*coefin/coefout))/coefin);
tabval=[cmin-0.0001:pas*10:cmax+0.0001];
tabval3=[cmin-0.0001:pas:cmax+0.0001];

cmap=jet(length(tabval3)-1); 
%  il=length(tabval)-1;
%  il5=round(il*8/100)
%  switch paramname
%      case 'PSAL'
%        cmap=cmocean('haline',il5+length(tabval)-1);  
%        cmap=cmap(il5+1:end,:);
%      case {'TEMP' , 'TPOT'} 
%          cmap=cmocean('thermal',il5+length(tabval)-1);  
%          cmap=cmap(il5+1:end,:);
%      case 'SIG0'
%           cmap=cmocean('dense',il5+length(tabval)-1);  
%          cmap=cmap(1:end-il5-1,:);
%  end 
hf=figure;
hold on
colormap(cmap)
set(gca,'fontsize',18)

switch vshading
    case 'faceted'
        pcolor(ones(nz+1,1)*double([cycnum;cycnum(end)+1]'),tabpres,tabparam)
        shading faceted
        set(gca,'xlim',[cycnum(1) cycnum(end)+1])
        xt=get(gca,'xtick');
        xtl=get(gca,'xticklabel');
        set(gca,'xtick',xt+0.5,'xticklabel',xtl)
    case 'interp'
%         set(gca,'xlim',[cycnum(1) cycnum(end)])
%         [c,h]=contourf(ones(nz+1,1)*double([cycnum;cycnum(end)]'),tabpres,tabparam,tabval);
         %keyboard
        pcolor(ones(nz+1,1)*double([cycnum;cycnum(end)+1]'),tabpres,tabparam)
        %keyboard
        %[c,h]=contourf(ones(nz+1,1)*double([cycnum;cycnum(end)]'),tabpres,tabparam,tabval);
        [c,h]=contour(ones(nz+1,1)*double([cycnum;cycnum(end)]'),tabpres,tabparam,tabval2,'k','LineStyle',':');
        [c,h]=contour(ones(nz+1,1)*double([cycnum;cycnum(end)]'),tabpres,tabparam,tabval,'k');
        shading flat
        set(gca,'xlim',[cycnum(1) cycnum(end)+1])
        xt=get(gca,'xtick');
        xtl=get(gca,'xticklabel');
        set(gca,'xtick',xt+0.5,'xticklabel',xtl)
        if isempty(strfind(paramname,'SIG0'))
        minsig=floor(min(min(tabsig(tabpres>500)))*10)/10;
        maxsig=ceil(max(max(tabsig(tabpres>500)))*10)/10;
        passig=(maxsig-minsig)/10;
        [c,h]= contour(ones(nz+1,1)*double([cycnum;cycnum(end)]'),tabpres,tabsig,[minsig:passig:maxsig],'w','LineWidth',1);
        clabel(c,h,'LabelSpacing',144*4,'FontSize',8,'Color','w');
        end
   
    case 'flat'
        pcolor(ones(nz+1,1)*double([cycnum;cycnum(end)+1]'),tabpres,tabparam)
        [c,h]=contour(ones(nz+1,1)*double([cycnum;cycnum(end)]'),tabpres,tabparam,tabval2,'Color',[0.4 0.4 0.4],'LineWidth',0.5);
        [c,h]=contour(ones(nz+1,1)*double([cycnum;cycnum(end)]'),tabpres,tabparam,tabval,'k','LineWidth',1);

        shading flat
        set(gca,'xlim',[cycnum(1) cycnum(end)+1])
        xt=get(gca,'xtick');
        xtl=get(gca,'xticklabel');
        set(gca,'xtick',xt+0.5,'xticklabel',xtl)
end
colorbar
caxis([cmin cmax]);
set(gca,'ydir','reverse','ylim',[200 10*max(ceil(pres(:)/10))])
%keyboard
box on
end


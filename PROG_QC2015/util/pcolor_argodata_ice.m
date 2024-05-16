function pcolor_argodata(cycnum,pres,param,paramname,vshading);
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
nprf=max(cycnum);
% interpolation sur une grille pres reguliere.
regpress=repmat([1:1:max(max(pres))]',1,nprf);

for ik=1:nprf
    icy=find(cycnum==ik);
    if isempty(icy)==0
        [a,b]=unique(pres(:,icy));   % correction cc 29/11/2021 : bug si pres n'est pas unique
        presuniq=unique(pres(:,icy));
        
        paramuniq=param(:,icy);
        iin=~isnan(presuniq);
        if sum(iin)>=2
            regparam(:,ik)=interp1(presuniq(iin),paramuniq(iin),regpress(:,ik));
        else
            regparam(:,ik)=NaN*regpress(:,ik);
        end
    else
        regparam(:,ik)=NaN*regpress(:,ik);
    end
end
param=regparam;
pres=regpress;
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

    


% defintion des intervalles
switch paramname
    case 'PSAL'
       coefin=10;
       coefout=1;
       pas=0.05;
    case {'TEMP' , 'TPOT'} 
        coefin=10;
        coefout=1;
        pas=0.1;
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
        pas=0.1;    
    case 'PSAT'
        coefin=1;
        coefout=10;
        pas=5;
end
cmin=coefout*min(floor(param(:)*coefin/coefout))/coefin;
cmax=coefout*max(ceil(param(:)*coefin/coefout))/coefin;
tabval=[cmin-0.001:pas:cmax+0.001];
cmap=jet(length(tabval)-1); 
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
%hf=figure;
hold on
box on
grid on
colormap(cmap)
set(gca,'fontsize',10)
%pcolor(ones(nz+1,1)*double([cycnum;cycnum(end)+1]'),tabpres,tabparam)
pcolor(ones(nz+1,1)*double([1:nprf+1]),tabpres,tabparam)
switch vshading
    case 'faceted'
        shading faceted
        set(gca,'xlim',[cycnum(1) cycnum(end)+1])
        xt=get(gca,'xtick');
        xtl=get(gca,'xticklabel');
        set(gca,'xtick',xt+0.5,'xticklabel',xtl)
    case 'interp'
        shading flat
        set(gca,'xlim',[cycnum(1) cycnum(end)])
        [c,h]=contour(ones(nz+1,1)*double([1:nprf+1]),tabpres,tabparam,tabval,'k');
    case 'flat'
        shading flat
        set(gca,'xlim',[cycnum(1) cycnum(end)+1])
        xt=get(gca,'xtick');
        xtl=get(gca,'xticklabel');
        set(gca,'xtick',xt+0.5,'xticklabel',xtl)
        %contourf(ones(nz+1,1)*double([1:nprf+1]),tabpres,tabparam,tabval);
        [c,h]=contour(ones(nz+1,1)*double([1:nprf+1]),tabpres,tabparam,[-1.8 -1.8],'w');
        clabel(c,h,'LabelSpacing',144*4,'FontSize',7,'Color','w','FontWeight','bold');
end
colorbar('horiz')
caxis([cmin cmax]);
set(gca,'ydir','reverse','ylim',[0 10*max(ceil(pres(:)/10))])


end


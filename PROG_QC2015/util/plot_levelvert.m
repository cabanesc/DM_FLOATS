function plot_levelvert(plotpath,wmonum,cycnum,pres,presqc,titflag,titsavedata)


[ncyc,nlev]=size(pres);

% détermination de newpres
% --------------------------------------------------------------
% Virginie THIERRY - 11/01/2008
% all pressure measurements at the same index levels

tabnpt=NaN*ones(size(cycnum,1),1);

for J = 1:size(cycnum,1)
    presvec   = pres(J,find(isfinite(pres(J,:))));
    indlev    = find(pres(J,:)==presvec(end)) ;
    tabnpt(J) = indlev(1);
end


[indlevmax,cyclevmax]=max(tabnpt);
vecpresref=pres(cyclevmax,:);
isf=find(isfinite(vecpresref));
vecpresref=vecpresref(isf);

newpres       = NaN*ones(size(cycnum,1),indlevmax);
newpresflag  = NaN*ones(size(cycnum,1),indlevmax);

presref=ones(size(cycnum,1),1)*vecpresref;


% presref est dimensionne a pres - valeurs flaguees a 4 pour 2047
% donc pas meme taille que newpres ...

for J=1:size(cycnum,1)
    
    if J == cyclevmax

        newpres(J,:) = pres(J,1:indlevmax);
        newpresflag(J,:)=presqc(J,1:indlevmax);
        presref(J,:) = pres(J,1:indlevmax);
       
    else
        presvec     = pres(J,find(isfinite(pres(J,:))));
        presvecflag = presqc(J,find(isfinite(pres(J,:))));
       
        npres=length(presvec);
        
        for jk=1:npres
            absdiff=abs(vecpresref-presvec(jk));
            vmin=min(absdiff);
            [indmin]=find(absdiff==vmin);
            if length(indmin)==2
                indlev=indmin(2);
            else
                indlev=indmin;
            end
            
             newpres(J,indlev)     = presvec(jk);
             newpresflag(J,indlev)= presvecflag(jk);
             presref(J,indlev)     = presvec(jk);
             
        end
    end
end

clear indlevmax cyclevmax tabnpt indlev

keyboard
tabval=newpres;
[~,nlev2] = size(tabval);
pas=50;
nomval='PRES';

vmin=0;
vmax=2100;
tabcont=[vmin:pas:vmax];
figure
orient landscape
set(gca,'fontsize',18)
colormap(jet(length(tabcont) -1 ));
cycnum=double(cycnum);
pcolor(cycnum*ones(1,nlev2),ones(ncyc,1)*[1:nlev2],tabval)
colorbar
caxis([vmin vmax])
set(gca,'ydir','reverse')
ylabel('Vertical level index')
xlabel('Cycle number')
title(['Float WMO ' wmonum ' - ' nomval]);
eval(['print -depsc2 ' plotpath wmonum '_' nomval '_raw_' titflag titsavedata '.eps']);

   
function plot_dbrut(plotpath,wmonum,paramx,paramy,nomparamx,nomparamy,labelx,labely,cycnum,lon,titflag)

% trace des donnees brutes en fonction de la pression
% param possibles : TEMP, PSAL, DOXY, SIG0, TPOT 
% cas particulier : PSAT


if  strcmp(nomparamx,'O2sat') == 0

    map=jet(length(cycnum));
 
    figure; 
    set(gca,'fontsize',18)
    grid
    hold on
    
    vymin = 0;
    vymax= 10*max(ceil(paramy(:)/10));
 
    for ii=1:length(lon)
        h=plot(paramx(ii,:),paramy(ii,:),'-');
        set(h,'color',map(ii,:));
        set(gca,'ydir','reverse');
        set(gca,'ylim',[vymin vymax]);

 %       set(gca,'ydir','reverse') 
    end
    
    title(['WMO ' wmonum ' - ' nomparamx],'interpreter','none');
    xlabel(labelx)
    ylabel(labely)
    eval(['print -depsc2 ' plotpath wmonum '_' nomparamy '_' nomparamx '_' titflag '.eps']);
    eval(['print -dpng ' plotpath wmonum '_' nomparamy '_' nomparamx '_' titflag '.png']);

else
    
    
% traces PSAT


    if length(cycnum)>2
        % Plot saturation aux 3 premiers niveaux

        cycnum2=cycnum;
        icz=find(cycnum==0);
        if length(icz) == 2
            cycnum2(icz(1))=-1;
        end

        map=jet(4);
        figure
        subplot(2,1,1),
        set(gca,'fontsize',18)
        grid
        hold on
 % revoir
        psat=paramx';
        for ii=1:4
            h=plot(cycnum2,psat(ii,:),'k-*');
            set(h,'color',map(ii,:))
        end
        ylabel(labelx);
        xlabel('Cycle number');
        set(gca,'xlim',[-1 max(cycnum2)+1]);
        subplot(2,1,2),
        set(gca,'fontsize',18)
        grid
        hold on
        % revoir
        pres = paramy';
        for ii=1:4
            h=plot(cycnum2,pres(ii,:),'k-*');
            set(gca,'ydir','reverse')
            set(h,'color',map(ii,:))
        end
        ylabel(labely);
        xlabel('Cycle number');
        set(gca,'xlim',[-1 max(cycnum2)+1]);

        hs=suptitle(['WMO ' wmonum ]);
        set(hs,'fontsize',18)
        eval(['print -depsc2 ' plotpath '/' wmonum '_' nomparamx '_' titflag '.eps']);
        eval(['print -dpng ' plotpath '/' wmonum '_' nomparamx '_' titflag '.png']);
    end
end
    


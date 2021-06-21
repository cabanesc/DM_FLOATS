function plot_inverdens(CONFIG,wmonum,cycnum,pres,temp,psal,titflag)


plotpath = [CONFIG.DIR_PLOT '/density_anomaly/'];
if ~exist(plotpath,'dir')
    mkdir(plotpath);
end


texte1= [' - Profil(s) avec une inversion de densite > 0.03  : '] ;

[nprf,nz]=size(pres);
float = num2str(wmonum);

dirsauv = [plotpath float];
if ~exist(dirsauv)
    eval(['mkdir ' dirsauv]);
end

filename=[dirsauv '/' float '_chkinv_sig' titflag '.txt'];
sauv = fopen(filename, 'w');

fprintf(sauv,'%s\n',texte1);
disp (texte1);

texte = ['Float ' float ];
fprintf(sauv,'%s\n',texte);
disp (texte)
texte= [datestr(now) ' - profile ' num2str(cycnum(1)) ' to ' num2str(cycnum(end))];
fprintf(sauv,'%s\n',texte);
disp (texte)


for icyc=1:nprf
    %keyboard
    [delta_upbot,delta_botup,is_inv_dens,rhop_i,rhop_ip,Pref] = test_densite14(pres(icyc,:),temp(icyc,:),psal(icyc,:)) ;
    
    niv_inv=find(is_inv_dens == 1);
    for ip=1:length(niv_inv)
        texte_sauve = [   'Cycle ' num2str(cycnum(icyc)) ' - Niv ' num2str(niv_inv(ip)) ' (P=' num2str(pres(icyc,niv_inv(ip))) ')'];
        fprintf(sauv,'%s\n',texte_sauve);
        disp (texte_sauve)
        texte_sauve = [   '..... upbot ' num2str(delta_upbot(niv_inv(ip))) ' - botup ' num2str(delta_botup(niv_inv(ip)))];
        fprintf(sauv,'%s\n',texte_sauve);
        disp(texte_sauve)
        
        
        h=figure;
        hold on
        title(['WMO ' float ' - Cycle ' num2str(cycnum(icyc)) ' - Niveau ' num2str(niv_inv(ip))  ' - Ecart ' num2str(rhop_i(niv_inv(ip))-rhop_ip(niv_inv(ip)))]);
        ideb = niv_inv(ip)-10;
        ifin = niv_inv(ip)+10;
        
        if  niv_inv(ip) < 21
            ideb = 1;
        end
        if  niv_inv(ip) > length(Pref)-10
            ifin = length(Pref);
        end
        grid
        plot(rhop_i(ideb:ifin),-pres(icyc,ideb:ifin),'.g',rhop_i(niv_inv(ip)),-pres(icyc,niv_inv(ip)),'.r','MarkerSize',17);
        plot(rhop_ip(ideb:ifin),-pres(icyc,ideb+1:ifin+1),'.c',rhop_ip(niv_inv(ip)),-pres(icyc,niv_inv(ip)+1),'*r','MarkerSize',17);
        xlabel('Potential density');
        ylabel('Pressure (db)');
        
        eval(['print -dpng ' plotpath float '/' float '_cyc' num2str(cycnum(icyc)) '_niv' num2str(niv_inv(ip)) '_densityanom.png']);
        
        close(h);
    end
end

end



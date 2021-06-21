function [bad_point]=select_interac_on_map(h,PLOT,bad_point)


ibad=size(bad_point,1);
dcm_obj = datacursormode(h); 
set(dcm_obj,'DisplayStyle','datatip',...
'SnapToDataVertex','off','Enable','on')
answer2=1;
answer3=menu(['Voulez vous afficher les points sucpects deja notés'],'YES','NO');
if answer3==1
all_bad_wmo = unique(bad_point(:,1));
for k=1:length(all_bad_wmo)
    iik=bad_point(:,1)==all_bad_wmo(k);
    iprofbad=bad_point(iik,2);
    jbad=(PLOT.wmobox.data==all_bad_wmo(k))&ismember(PLOT.profil.data,iprofbad);
    %scatter(PLOT.longitude.data(jbad),PLOT.latitude.data(jbad),40,PLOT.psal.data(jbad),'filled');
    scatter(PLOT.longitude.data(jbad),PLOT.latitude.data(jbad),40,'m','filled');
end
end
while (answer2==1|answer2==2)
    answer2=menu(['Click on a point on the figure and then  ADD or REMOVE it from the list of bad points'],'ADD','REMOVE','STOP');
    cursor_info = getCursorInfo(dcm_obj);
    ipos=find(cursor_info.Position(1)==PLOT.longitude.data&cursor_info.Position(2)==PLOT.latitude.data);
    
    [thebad_prof]=PLOT.profil.data(ipos)
    thebad_wmo=PLOT.wmobox.data(ipos)
    PLOT.psal.data(ipos)
    for kk=1:length(thebad_prof)
        if answer2==1
            ibad = ibad+1;
            bad_point(ibad,1)=thebad_wmo(kk);
            bad_point(ibad,2)=thebad_prof(kk);
            figure(h)
            scatter(PLOT.longitude.data(ipos),PLOT.latitude.data(ipos),40,'m','filled');
            
        elseif answer2==2
            jj = find((bad_point(:,1)==thebad_wmo(kk) & bad_point(:,2)==thebad_prof(kk)));
            if isempty(jj)==0
                ibad=ibad-length(jj);
                bad_point(jj,:)=[];
            end
            figure(h)
            scatter(PLOT.longitude.data(ipos),PLOT.latitude.data(ipos),40,PLOT.psal.data(ipos),'filled');
        end
    end

end
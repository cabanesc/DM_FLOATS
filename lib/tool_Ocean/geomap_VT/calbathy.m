load /home1/capchat/reference/AN_smithsandwell
[nlat nlon]=size(bathy);
indj=0;
for ij=6:5:nlat-5
  indj=indj+1;
  lat_bathy2(indj)=mean(lat_bathy([ij-5:ij+5]));
end
indi=0;
for ii=6:5:nlon-5
  indi=indi+1;
  lon_bathy2(indi)=mean(lon_bathy([ii-5:ii+5]));
end

indi=0;
for ii=6:5:nlon-5
  indi=indi+1;
  indj=0;
  for ij=6:5:nlat-5
    indj=indj+1;
    tmp=bathy(ij-5:ij+5,ii-5:ii+5);
    bathy2(indj,indi)=mean(tmp(:));
  end
end
[nlat nlon]=size(bathy2);

figure
m_proj('lambert','long',[-45 0],...
                 'lat',[35 65]);
%m_gshhs_i('color','k');
m_grid('box','fancy','tickdir','in');	     
hold on

v=[-2000 -2000];   
[n,p]=m_contour(lon_bathy2,lat_bathy2,bathy2,v,'color',[.6 .6 .6]);
v=[-1000 -1000];
[q,r]=m_contour(lon_bathy2,lat_bathy2,bathy2,v,'color',[.7 .7 .7]);
v=[-200 -200]; 
[c,h]=m_contour(lon_bathy2,lat_bathy2,bathy2,v,'color',[.8 .8 .8]);
v=[0 0]; 
[c,h]=m_contour(lon_bathy2,lat_bathy2,bathy2,v,'k','linewidth',2);


save AN_smithsandwell_lowres bathy2 lat_bathy2 lon_bathy2
function [h1,h2,h3]=gebco(cart,fichier,zone,bathy,type_cont,type_map,beta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Trace la bathymetrie a partir d'un fichier GEBCO
% avec une projection de Mercator sur une ellipsoide WGS-84
%
% V1.0 Y.A. 02/01/96 : Creation
%
% function [h1,h2,h3]=gebco(cart,fichier,zone,bathy,type_cont,type_map,beta)
%
% Si cart=0: Trace de la carte 
%           fichier  = nom du fichier ou des fichiers GEBCO
%           zone     = zone de la carte ex = [28.0,38.0,-26.0,-16.0]
%                      (valeurs entieres)
%           bathy    = tableau contenant les bathy par ordre croissant
%                      ex: [0,500,1000,1500,2000,2500]
%           type_cont= type de contour.
%                      0 = Fill (a utiliser si les bathys sont fermees)
%                      1 = Line (a utiliser si les bathys sont ouvertes)
%           type_map = Type de colormap
%                      1  = bone
%                      2  = cool
%                      3  = copper
%                      4  = flag
%                      5  = gray
%                      6  = hot
%                      7  = hsv
%                      8  = jet
%                      9  = pink
%                      10 = prism
%           beta     = Luminosite
%                      0  normal
%                      0   < beta <= 1 plus clair
%                      -1 <= beta <  0 plus clair
%
% Si cart=1: Trace du cartouche 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if (cart==0)

%
% Dessin de la fenetre
%
ratio = (zone(4)-zone(3))/(proj_mercator(zone(2))-proj_mercator(zone(1)));

h=figure('Units','normal','Position',[0.1,0.1,0.8,0.8]);
h1=[h1,h];
hold on;
axis('off');

h=axes('Position',[0.1,0.1,0.4,0.8],'Visible','off','AspectRatio',[ratio,1],'Xlim',[zone(3),zone(4)],'Ylim',[proj_mercator(zone(1)),proj_mercator(zone(2))]);
h1=[h1,h];
hold on;

%
% Calcul de la palette de couleur
%
if (type_cont==0)
	niv_coul=length(bathy)+1;
else
	niv_coul=length(bathy);
end
if (type_map==1)  map=bone(niv_coul);   end
if (type_map==2)  map=cool(niv_coul);   end
if (type_map==3)  map=copper(niv_coul); end
if (type_map==4)  map=flag(niv_coul);   end
if (type_map==5)  map=gray(niv_coul);   end
if (type_map==6)  map=hot(niv_coul);    end
if (type_map==7)  map=hsv(niv_coul);    end
if (type_map==8)  map=jet(niv_coul);    end
if (type_map==9)  map=pink(niv_coul);   end
if (type_map==10) map=prism(niv_coul);  end
map(niv_coul,1:3)=[.5 .5 .5];
colormap(map);
if (beta~=0) brighten(beta); end
map=colormap;
hold on;

%
% Remplissage du fond de carte
%
if (type_cont==0)
	h=fill([zone(3),zone(3),zone(4),zone(4)],[proj_mercator(zone(1)),proj_mercator(zone(2)),proj_mercator(zone(2)),proj_mercator(zone(1))],map(1,:));
	hold on;
end

%
% contourage de la bathy (Mode fill)
%
if (type_cont==0)
	for index_fic=1:length(fichier(:,1))
		for index=length(bathy):-1:1
			flag=0;
			fid=fopen(fichier(index_fic,:));
			ente=fscanf(fid,'%g %g',[2 1]);
			trac=fscanf(fid,'%g %g',[2 ente(2)]);
			while(ente(1)<bathy(index))
				ente=fscanf(fid,'%g %g',[2 1]);
				if (length(ente)>0) trac=fscanf(fid,'%g %g',[2 ente(2)]); else flag=1; break; end
			end
			if (flag==0)
				while(ente(1)==bathy(index))
					h=fill(trac(2,:),proj_mercator(trac(1,:)),map(niv_coul-index+1,:));
					hold on;
					ente=fscanf(fid,'%g %g',[2 1]);
					if (length(ente)>0) trac=fscanf(fid,'%g %g',[2 ente(2)]); else break; end
				end
			end
			fclose(fid);
		end
	end
else
	for index_fic=1:length(fichier(:,1))
		flag=0;
		fid=fopen(fichier(index_fic,:));
		ente=fscanf(fid,'%g %g',[2 1]);
		trac=fscanf(fid,'%g %g',[2 ente(2)]);
		for index=1:length(bathy)
			if (flag==0)
				while(ente(1)<bathy(index))
					ente=fscanf(fid,'%g %g',[2 1]);
					if (length(ente)>0) trac=fscanf(fid,'%g %g',[2 ente(2)]); else flag=1; break; end
				end
			end
			if (flag==0)
				while(ente(1)==bathy(index))
					h=line(trac(2,:),proj_mercator(trac(1,:)),'Color',map(niv_coul-index+1,:));
					hold on;
					ente=fscanf(fid,'%g %g',[2 1]);
					if (length(ente)>0) trac=fscanf(fid,'%g %g',[2 ente(2)]); else flag=1; break; end
				end
			end
		end
		fclose(fid);
	end
end
orient landscape;

%
% Trace de la grille et des graduations
%
for index=zone(1):1:zone(2)
	trx=[zone(3), zone(4)];
	tr_y=[proj_mercator(index), proj_mercator(index)];
	if (type_cont==0)
		h=plot(trx,tr_y,'k:','linewidth',1);
	else
		h=plot(trx,tr_y,'w:','linewidth',1);
	end
	hold on;
	texte='0';
	if (index > 0) texte=['N',num2str(abs(index))]; end
	if (index < 0) texte=['S',num2str(abs(index))]; end
	h=text(zone(3)-.2,proj_mercator(index),texte,'fontsize',10,'Vert','middle','Hor','right');
	h=text(zone(4)+.2,proj_mercator(index),texte,'fontsize',10,'Vert','middle','Hor','left');
	hold on;
end
for index=zone(3):1:zone(4)
	tr_y=[proj_mercator(zone(1)), proj_mercator(zone(2))];
	trx=[index, index];
	if (type_cont==0)
		h=plot(trx,tr_y,'k:','linewidth',1);
	else
		h=plot(trx,tr_y,'w:','linewidth',1);
	end
	hold on;
	texte='0';
	if (index > 0) texte=['E',num2str(abs(index))]; end
	if (index < 0) texte=['W',num2str(abs(index))]; end
	h=text(index,proj_mercator(zone(1))-.2,texte,'fontsize',10,'Vert','top','Hor','center');
	h=text(index,proj_mercator(zone(2))+.2,texte,'fontsize',10,'Vert','bottom','Hor','center');
	hold on;
end
for index=1:2 
	trx=zone(3):zone(4);
	tr_y=proj_mercator(zone(index))*ones(max(size(trx)),1)';
	h=plot(trx,tr_y,'w','linewidth',2);
	hold on;
end
for index=3:4
	tr_y=[proj_mercator(zone(1)), proj_mercator(zone(2))];
	trx=[zone(index), zone(index)];
	h=plot(trx,tr_y,'w','linewidth',2);
	hold on;
end

save gebco_tmp.mat;

return
end


%
% Construction de la ColorBar
%

load gebco_tmp.mat;

h=axes('Position',[0.56,0.1,0.02,0.8],'Visible','off','Xlim',[0,1],'Ylim',[0,1]);
h2=[h2,h];
hold on;

for index = niv_coul:-1:1
	x=[0 1];
	y=[((index/niv_coul)-1/niv_coul) index/niv_coul];
	h=image(x,y,index);
	hold on;
	if (type_cont==0)
		if (index==niv_coul)
			if (bathy(1)~=0)
				texte=['supérieur à ' num2str(-bathy(1)) 'm'];
			else
				texte=['supérieur à 0m'];
			end
		end
		if (index>1)&(index<niv_coul) 
			if (bathy(niv_coul-index)~=0)
				texte=[num2str(-bathy(niv_coul-index+1)) 'm à ',num2str(-bathy(niv_coul-index)) 'm'];
			else
				texte=[num2str(-bathy(niv_coul-index+1)) 'm à 0m'];
			end
		end
		if (index==1) texte=['inférieur à ' num2str(-bathy(niv_coul-1)) 'm']; end
	else
		if (bathy(niv_coul-index+1)~=0)
			texte=[num2str(-bathy(niv_coul-index+1)) 'm'];
		else
			texte=['0m'];
		end
	end
	h=text(1.2,((index/niv_coul)-0.5/niv_coul),texte,'fontsize',8,'Vert','middle','Hor','left');
	hold on;
end

trx=[0,1];
tr_y=[0,0];
h=plot(trx,tr_y,'w-','linewidth',1);
hold on;

trx=[0,1];
tr_y=[1,1];
h=plot(trx,tr_y,'w-','linewidth',1);
hold on;

trx=[0,0];
tr_y=[0,1];
h=plot(trx,tr_y,'w-','linewidth',1);
hold on;

trx=[1,1];
tr_y=[0,1];
h=plot(trx,tr_y,'w-','linewidth',1);
hold on;


%
% Construction du cartouche
%
h=axes('Position',[0.7,0.1,0.2,0.8],'Visible','off','Xlim',[0,1],'Ylim',[0,1]);
h3=[h3,h];
hold on;

trx=[0,1];
tr_y=[0,0];
h=plot(trx,tr_y,'w-','linewidth',1);
hold on;

trx=[0,1];
tr_y=[0.5,0.5];
h=plot(trx,tr_y,'w-','linewidth',1);
hold on;

trx=[0,1];
tr_y=[1,1];
h=plot(trx,tr_y,'w-','linewidth',1);
hold on;

trx=[0,0];
tr_y=[0,1];
h=plot(trx,tr_y,'w-','linewidth',1);
hold on;

trx=[1,1];
tr_y=[0,1];
h=plot(trx,tr_y,'w-','linewidth',1);
hold on;

h=text(0.5,0.375,['Projection: MERCATOR'],'fontsize',10,'Vert','middle','Hor','center');
h=text(0.5,0.250,['Ellipsoide: WGS-84'],'fontsize',10,'Vert','middle','Hor','center');
h=text(0.5,0.125,['IFREMER - LPO'],'fontsize',12,'Vert','middle','Hor','center');

delete gebco_tmp.mat

return


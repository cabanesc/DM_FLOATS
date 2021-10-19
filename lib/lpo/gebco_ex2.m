%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Trace de la zone Thetis2   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Initialisation
%clear;
%close all;
%
%% fichiers GEBCO de la zone
%fichier(1,:)='/home/corduan/yauffret/ibcm/ibcm01.asc';
%fichier(2,:)='/home/corduan/yauffret/ibcm/ibcm02.asc';
%fichier(3,:)='/home/corduan/yauffret/ibcm/ibcm06.asc';
%fichier(4,:)='/home/corduan/yauffret/ibcm/ibcm07.asc';
%
%% trace de la zone 36N-45N et 1E-9E avec la bathy de 0 a 2800m tous les 200m
%% la palette utilisee est la 'jet' avec une correction de luminosite de 0.8
%gebco(0,fichier,[36.0,45.0,1.0,9.0],[0:200:2800],1,8,0.7);
%
%% Ajout des instruments dans le trace
%W1=[proj_mercator(39.57),1.72];
%W2=[proj_mercator(37.24),4.29];
%W3=[proj_mercator(37.35),6.80];
%W4=[proj_mercator(40.50),7.45];
%W5=[proj_mercator(42.52),8.37];
%H =[proj_mercator(42.65),5.19];
%S =[proj_mercator(39.50),4.75];
%
%plot([W1(2) W2(2)],[W1(1) W2(1)],'r-');
%plot([W2(2) W3(2)],[W2(1) W3(1)],'r-');
%plot([W3(2)  H(2)],[W3(1)  H(1)],'r-');
%plot([W2(2)  H(2)],[W2(1)  H(1)],'r-');
%plot([W1(2)  H(2)],[W1(1)  H(1)],'r-');
%
%plot([W5(2) W1(2)],[W5(1) W1(1)],'g-.');
%plot([W5(2) W2(2)],[W5(1) W2(1)],'g-.');
%plot([W5(2) W3(2)],[W5(1) W3(1)],'g-.');
%plot([W4(2) W2(2)],[W4(1) W2(1)],'g-.');
%
%plot([H(2)  W5(2)],[H(1)  W5(1)],'y-');
%plot([W5(2)  S(2)],[W5(1)  S(1)],'y-');
%plot([W3(2)  S(2)],[W3(1)  S(1)],'y-');
%plot([W2(2)  S(2)],[W2(1)  S(1)],'y-');
%
%plot([H(2)   S(2)],[H(1)   S(1)],'m--');
%plot([H(2)  W4(2)],[H(1)  W4(1)],'m--');
%plot([S(2)  W4(2)],[S(1)  W4(1)],'m--');
%plot([W3(2) W4(2)],[W3(1) W4(1)],'m--');
%plot([W5(2) W4(2)],[W5(1) W4(1)],'m--');
%
%plot(W1(2),W1(1),'w*');
%text(W1(2),W1(1),'W1','Color','w');
%plot(W2(2),W2(1),'w*');
%text(W2(2),W2(1),'W2','Color','w');
%plot(W3(2),W3(1),'w*');
%text(W3(2),W3(1),'W3','Color','w');
%plot(W4(2),W4(1),'w*');
%text(W4(2),W4(1),'W4','Color','w');
%plot(W5(2),W5(1),'w*');
%text(W5(2),W5(1),'W5','Color','w');
%plot(S(2) ,S(1) ,'w*');
%text(S(2) ,S(1) ,'S','Color','w');
%plot(H(2) ,H(1) ,'w*');
%text(H(2) ,H(1) ,'H','Color','w');
%
%% trace du cartouche
%gebco(1);
%
%% Ajout de commentaires et legende  dans le cartouche
%text(0.5,0.95,['THETIS 2'],'fontsize',16,'Vert','middle','Hor','center');
%
%text(0.5,0.90,['Données IBCM'],'fontsize',6,'Vert','middle','Hor','center');
%
%plot([0.1 0.35],[0.85 0.85],'y-');
%text(0.35,0.85,[': réciproque'],'fontsize',10,'Vert','middle','Hor','left');
%
%plot([0.1 0.35],[0.75 0.75],'r-');
%text(0.35,0.75,[': %semi-réciproque'],'fontsize',10,'Vert','middle','Hor','left');
%
%plot([0.1 0.35],[0.65 0.65],'m--');
%text(0.35,0.65,[': simple'],'fontsize',10,'Vert','middle','Hor','left');
%
%plot([0.1 0.35],[0.55 0.55],'g-.');
%text(0.35,0.55,[': trace simple'],'fontsize',10,'Vert','middle','Hor','left');


% Initialisation
clear;
close all;

% fichiers GEBCO de la zone
fichier(1,:)='/home/corduan/yauffret/ibcm/ibcm01.asc';
fichier(2,:)='/home/corduan/yauffret/ibcm/ibcm02.asc';
fichier(3,:)='/home/corduan/yauffret/ibcm/ibcm06.asc';
fichier(4,:)='/home/corduan/yauffret/ibcm/ibcm07.asc';

% trace de la zone 36N-45N et 1E-9E avec la bathy de 0 a 2800m tous les 200m
% la palette utilisee est la 'jet' avec une correction de luminosite de 0.8
gebco(0,fichier,[36.0,45.0,1.0,9.0],[0:200:2800],1,8,0.7);

% Ajout des instruments dans le trace
W1=[proj_mercator(39.57),1.72];
W2=[proj_mercator(37.24),4.29];
W3=[proj_mercator(37.35),6.80];
W4=[proj_mercator(40.50),7.45];
W5=[proj_mercator(42.52),8.37];
H =[proj_mercator(42.65),5.19];
S =[proj_mercator(39.50),4.75];

plot([W1(2) W2(2)],[W1(1) W2(1)],'r-');
plot([W2(2) W3(2)],[W2(1) W3(1)],'r-');
plot([W3(2)  H(2)],[W3(1)  H(1)],'r-');
plot([W2(2)  H(2)],[W2(1)  H(1)],'r-');
plot([W1(2)  H(2)],[W1(1)  H(1)],'r-');

plot([W5(2) W1(2)],[W5(1) W1(1)],'g-.');
plot([W5(2) W2(2)],[W5(1) W2(1)],'g-.');
plot([W5(2) W3(2)],[W5(1) W3(1)],'g-.');
plot([W4(2) W2(2)],[W4(1) W2(1)],'g-.');

plot([H(2)  W5(2)],[H(1)  W5(1)],'y-');
plot([W5(2)  S(2)],[W5(1)  S(1)],'y-');
plot([W3(2)  S(2)],[W3(1)  S(1)],'y-');
plot([W2(2)  S(2)],[W2(1)  S(1)],'y-');

plot([H(2)   S(2)],[H(1)   S(1)],'m--');
plot([H(2)  W4(2)],[H(1)  W4(1)],'m--');
plot([S(2)  W4(2)],[S(1)  W4(1)],'m--');
plot([W3(2) W4(2)],[W3(1) W4(1)],'m--');
plot([W5(2) W4(2)],[W5(1) W4(1)],'m--');

plot(W1(2),W1(1),'w*');
text(W1(2),W1(1),'W1','Color','w');
plot(W2(2),W2(1),'w*');
text(W2(2),W2(1),'W2','Color','w');
plot(W3(2),W3(1),'w*');
text(W3(2),W3(1),'W3','Color','w');
plot(W4(2),W4(1),'w*');
text(W4(2),W4(1),'W4','Color','w');
plot(W5(2),W5(1),'w*');
text(W5(2),W5(1),'W5','Color','w');
plot(S(2) ,S(1) ,'w*');
text(S(2) ,S(1) ,'S','Color','w');
plot(H(2) ,H(1) ,'w*');
text(H(2) ,H(1) ,'H','Color','w');

% trace du cartouche
gebco(1);

% Ajout de commentaires et legende  dans le cartouche
text(0.5,0.95,['THETIS 2'],'fontsize',16,'Vert','middle','Hor','center');

text(0.5,0.90,['Données IBCM'],'fontsize',6,'Vert','middle','Hor','center');

plot([0.1 0.35],[0.85 0.85],'y-');
text(0.35,0.85,[': réciproque'],'fontsize',10,'Vert','middle','Hor','left');

plot([0.1 0.35],[0.75 0.75],'r-');
text(0.35,0.75,[': semi-réciproque'],'fontsize',10,'Vert','middle','Hor','left');

plot([0.1 0.35],[0.65 0.65],'m--');
text(0.35,0.65,[': simple'],'fontsize',10,'Vert','middle','Hor','left');

plot([0.1 0.35],[0.55 0.55],'g-.');
text(0.35,0.55,[': trace simple'],'fontsize',10,'Vert','middle','Hor','left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Trace de la zone cambios   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Initialisation
%clear;
%close all;
%
%% trace de la zone 28N-38N et 26W-16W avec la bathy de 0 a 4500m tous les 500m
%% la palette utilisee est la 'jet' avec une correction de luminosite de 0.8
%gebco(0,'/home/corduan/yauffret/gebco/gebco508.asc',[28.0,38.0,-26.0,-16.0],[0:500:4500],0,8,0.6);
%
%% trace du cartouche
%gebco(1);
%
%% Ajout d'un commentaire dans le cartouche
%text(0.5,0.95,['CAMBIOS'],'FontSize',16,'Vert','middle','Hor','center');


% Initialisation
clear;
close all;

% trace de la zone 28N-38N et 26W-16W avec la bathy de 0 a 4500m tous les 500m
% la palette utilisee est la 'jet' avec une correction de luminosite de 0.8
gebco(0,'/home/corduan/yauffret/gebco/gebco508.asc',[28.0,38.0,-26.0,-16.0],[0:500:4500],0,8,0.6);

% trace du cartouche
gebco(1);

% Ajout d'un commentaire dans le cartouche
text(0.5,0.95,['CAMBIOS'],'FontSize',16,'Vert','middle','Hor','center');


%
% function [fcorio,beta]=fcor(xlat)
%
% calcule le parametre de coriolis et beta en fonction de la latitude
%
% fcorio = 2*omega * sin(theta)
% beta   = 2*omega * cos(theta) /R
%
% input: xlat (degres)
% outout: fcor et beta

function [fcorio,beta]=fcor(xlat)

R=6371*1000;
theta=2*pi*xlat/360; % degres
omega=2*pi/(24*3600); % rad/s



fcorio  = 2*omega*sin(theta);
beta    = 2*omega*cos(theta)/R;
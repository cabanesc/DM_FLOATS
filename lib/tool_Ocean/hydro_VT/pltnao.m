
function [yearind,naoind]=getnao()

inpath='/home/revellata/vthierry/ATLNORD/NAO';

load([inpath 'naodjfmindex.1864-2004.asc']);

yearind=naodjfmindex(:,1);
naoind=naodjfmindex(:,2);

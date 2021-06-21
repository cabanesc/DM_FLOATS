    
    
%------------------------------------------------------
% C.Lagadec - septembre 2015
% Lecture de tous les monoprofifs d'un flotteur
% restitution dans structure FLm 
% a� l'aide du programme ecrit par Cecile 
% et mise dans variables appropriees pour traitement
% verif_flag et plotdata
%------------------------------------------------------

pltoxy = 'n';

lat     = FLm.latitude.data;
lon     = FLm.longitude.data;
juld    = FLm.juld.data;
juld_qc = FLm.juld_qc.data;
position_qc = FLm.position_qc.data;
cycnum      = FLm.cycle_number.data; 

if datatype == 2
adj='_adjusted';
else
adj='';
end

platform_number = FLm.platform_number.data;
direction = FLm.direction.data;
pres = FLm.(['pres' adj]).data';
fillvalue=FLm.(['pres' adj]).FillValue_;

pres(find(pres==fillvalue))=NaN;

[nlev,ncyc] = size(pres);


temp = FLm.(['temp' adj]).data';
psal = FLm.(['psal' adj]).data'; 
pres_qc = FLm.(['pres' adj '_qc']).data';
temp_qc = FLm.(['temp' adj '_qc']).data';
psal_qc = FLm.(['psal' adj '_qc']).data';
psal(find(psal==fillvalue))=NaN;
temp(find(temp==fillvalue))=NaN;

tpot = sw_ptmp(psal,temp,pres,0);
sig0 = sw_pden(psal,temp,pres,0) - 1000;


%-------------------------------------------------------

if strcmp(pltoxy,'o')
    doxy = FLm.doxy.data';
    doxy(find(doxy==fillvalue))=NaN;
    doxy_qc=FLm.doxy_qc.data';
end

% transformation des flags en num

for ic = 1:ncyc
    for in=1:nlev    
        if strcmp(pltoxy,'o')
            doxy_qc_num(in,ic) = str2double(doxy_qc(in,ic));
        end
        
        pres_qc_num(in,ic)=str2double(pres_qc(in,ic));
        psal_qc_num(in,ic)=str2double(psal_qc(in,ic));
        temp_qc_num(in,ic)=str2double(temp_qc(in,ic));
        sig0_qc_num(in,ic)=max(temp_qc_num(in,ic),psal_qc_num(in,ic));
        
    end
end

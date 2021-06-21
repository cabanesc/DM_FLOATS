    function [PLOT,CTD] = find_psal_on_z(CTD, P_MIN, P_MAX, PLOT)
    
    if nargin <=4
    PLOT=[];
    end
    
    % trouve les données appartenant a un intervalle de pression [P_MIN : P_MAX]
    itpot=find(CTD.pres.data>=P_MIN & CTD.pres.data<=P_MAX);% & repmat(CTD.dates.data,1,size(CTD.pres.data,2))>=datenum(2000,01,01));
    masque=1*(CTD.pres.data>=P_MIN & CTD.pres.data<=P_MAX );%& repmat(CTD.dates.data,1,size(CTD.pres.data,2))>=datenum(2000,01,01));
    masque(masque==0)=NaN;
    CTD.psal_mean.data=meanoutnan(CTD.psal.data.*masque,2);
    CTD.pres_mean.data=meanoutnan(CTD.pres.data.*masque,2);
    CTD.tpot_mean.data=meanoutnan(CTD.tpot.data.*masque,2);
    notnan = (~isnan(CTD.psal_mean.data));
    % sauvegarde des données du plot.
    
    if isempty(PLOT)==1
    ifin1=0;
    else
    ifin1=length(PLOT.latitude.data)
    end
    
    ideb1 = ifin1+1;
    ifin1 = ideb1+sum(notnan)-1;
    
    PLOT.latitude.data(ideb1:ifin1) = CTD.latitude.data(notnan);
    PLOT.longitude.data(ideb1:ifin1) = CTD.longitude.data(notnan);
    PLOT.psal.data(ideb1:ifin1) = CTD.psal_mean.data(notnan);
    if isfield(CTD,'dates')
    PLOT.dates.data(ideb1:ifin1) = CTD.dates.data(notnan);
    end
    if isfield(CTD,'profil')
    PLOT.profil.data(ideb1:ifin1) = CTD.profil.data(notnan);
    end
    if isfield(CTD,'wmobox')
    PLOT.wmobox.data(ideb1:ifin1) = CTD.boxes.data(notnan);
    end
    return
function read_history(Co,Dim,NcVar,fw1)
% -========================================================
%   USAGE :   print_history(filename) single profile file
%   PURPOSE : Print the history of the QC performed on an Argo profile
% -----------------------------------
%   INPUT :
%     filename   (string)  - Argo file name-
% -----------------------------------
%   OUTPUT : printed
%
% -----------------------------------
%   HISTORY  : created (2013) ccabanes
%            : modified (yyyy) byxxx
%   CALLED SUBROUTINES: read_netcdf_allthefile
% ========================================================

% read the netcdf file

%[Co,Dim] = read_netcdf_allthefile(filename,NcVar);

Co=check_FirstDimArray_is(Co,'N_HISTORY');
%keyboard

if Dim.n_prof.dimlength>1
    %  disp('Multi profile file')
end
ifl=0;
for n_prof=1:Dim.n_prof.dimlength % pour chaque cycle
    profile_number=n_prof;
    icy=0;
    isprint=0;
	
    prev_start_pres=10000;
    prev_stop_pres=-10000;
    min_pres=10000;
    max_pres=-10000;
    
    % cherche les dates de l'historique
    for k=1:Dim.n_history.dimlength
        thedate = (squeeze(Co.history_date.data(k,profile_number,:)))';
        thedatenum(k) = datenum(thedate,'yyyymmddHHMMSS');
    end
    
    [sortdate,isort]=sort(thedatenum);
    
    for k=[isort] % pour chaque date trouvee
        % HISTORY_DATE
        thedate = (squeeze(Co.history_date.data(k,profile_number,:)))';
        thedatevec = datevec(thedate,'yyyymmddHHMMSS');
        thedatestr = datestr(thedatevec,'dd/mm/yyyy');
        
        % HISTORY_SOFTWARE
        software=squeeze(Co.history_software.data(k,profile_number,:))';
        
        if strcmp(strtrim(software),'COAR')
            software= 'COAR (real time tests) ';
        end
        
        if strcmp(strtrim(software),'COOA')
            software= 'COOA (Coriolis Objective Analysis test) ';
        end
        
        if strcmp(strtrim(software),'SCOO')
            software= 'SCOOP ';
        end
        
        % HISTORY_SOFTWARE_RELEASE
        software_release=squeeze(Co.history_software_release.data(k,profile_number,:))';
        if strcmp(strtrim(software_release),'1.4')
            software= 'SCOOP2 ';
        end
       
        action=squeeze(Co.history_action.data(k,profile_number,:))';
        
        %  traduction of the action codes  (reference table 7, Argo user manual)
        if strcmp(strtrim(action),'QCP$')
            action=' THESE TESTS WERE PERFORMED ';
        end
        
        if strcmp(strtrim(action),'QCF$')
            action= ' THESE TESTS FAILED ';
        end
        
        if strcmp(strtrim(action),'QC')
            action= ' VISUAL QUALITY CONTROL ';
        end
        
        if strcmp(strtrim(action),'CF')
            action= ' FLAG WERE CHANGED ';
        end
        
        
        
        
        % HISTORY_PARAMETER
        parameter=squeeze(Co.history_parameter.data(k,profile_number,:))';
        previous_value=squeeze(Co.history_previous_value.data(k,profile_number,:))';
        start_pres=squeeze(Co.history_start_pres.data(k,profile_number,:))';
        stop_pres=squeeze(Co.history_stop_pres.data(k,profile_number,:))';
        
        if previous_value==Co.history_previous_value.FillValue_
            previous_value=[];
        end
        if start_pres==Co.history_start_pres.FillValue_
            start_pres=[];
        end
        if stop_pres==Co.history_stop_pres.FillValue_
            stop_pres=[];
        end
       %keyboard
        
        
        % PRINT 

        cycle = [sprintf('%3.3i',Co.cycle_number.data(profile_number,:)), strtrim(Co.direction.data(profile_number,:))];
        if isempty(previous_value)==0&strcmp(strtrim(software),'SCOOP2')
            
            % recherche le nouveau qc dans PARAM_QC
            iniv = find(Co.pres.data(profile_number,:)==start_pres);
            new_flag = Co.([strtrim(lower(parameter)) '_qc']).data(profile_number,iniv);
            platform_number = strtrim(Co.platform_number.data(profile_number,:));
            
            if str2num(new_flag)~=3
                
                
                if ifl==0 % compteur mis a zero a chaque nouveau flotteur
                    strto_write=[strtrim(Co.platform_number.data(profile_number,:)) ' & ' cycle '&' strtrim(parameter) '&' num2str(previous_value) '&' new_flag '&' num2str(start_pres) ' : ' num2str(stop_pres) '&' thedatestr '\\'];
                    prev_strto_write=strto_write
                    save_min_pres=start_pres;
                    save_max_pres=stop_pres;
                   
                else
                    if icy==0 % compteur mis a zero a chaque nouveau cycle
                        strto_write=[ ' & ' cycle '&' strtrim(parameter) '&' num2str(previous_value) '&' new_flag '&' num2str(start_pres) ' : ' num2str(stop_pres)  '&' thedatestr '\\'];
                        prev_strto_write=strto_write
                        save_min_pres=start_pres;
                        save_max_pres=stop_pres;
                       
                    else
                        
                        strto_write=[ ' & &' strtrim(parameter) '&' num2str(previous_value) '&' new_flag '&' num2str(start_pres) ' : ' num2str(stop_pres)  '&' thedatestr '\\'];
                        
                        % if strcmp(strtrim(parameter),strtrim(prev_parameter))&prev_previous_value==previous_value&strcmp(new_flag,prev_new_flag)
                            % % meme flotteur, cycles,parametre,flag : on condense
                            % min_pres=min(start_pres,save_min_pres)
                            % max_pres=max(stop_pres,save_max_pres) 
                            % % ATTETION********** bug ici si  save_min_pres ou save_max_pres est egale a la valeur d'un flag (1 2 3 ou 4)!!
                            % % CORRECTION du 11/03/2011 -> trouve les séparateurs '&' et change seulement les pressions.
                            % iseparateur=strfind(prev_strto_write,'&');
                            
                            % old_str_pres = prev_strto_write(iseparateur(5)+1:iseparateur(6)-1);
                            % new_str_pres = [num2str(start_pres) ' : ' num2str(stop_pres)];
                            % strto_write=strrep(prev_strto_write,old_str_pres,new_str_pres);
                            
                            % save_min_pres=min_pres;
                            % save_max_pres=max_pres;
                            % %strto_write=[ ' & &' strtrim(parameter) '&' num2str(previous_value) '&' new_flag '&' num2str(min_pres) ' : ' num2str(max_pres)  '&' thedatestr '\\'];
                            % %thesame=1;
                        %else
                            save_min_pres=start_pres;
                            save_max_pres=stop_pres;
                            fprintf(fw1,'%s\n',prev_strto_write)
                            isprint==isprint+1;;
                        %end
                        prev_strto_write=strto_write;
                        
                        
                        %fprintf(fw1,'%s\n', [ ' & &' strtrim(parameter) '&' num2str(previous_value) '&' new_flag '&' num2str(start_pres) ' : ' num2str(stop_pres)  '&' thedatestr '\\']);
                    end
                end
                icy=icy+1;
                ifl=ifl+1;
                prev_platform_number=Co.platform_number.data(profile_number,:);
                prev_iniv=iniv;
                prev_new_flag=new_flag;
                prev_cycle=cycle;
                prev_parameter=parameter;
                prev_previous_value=previous_value;
                prev_start_pres=start_pres;
                prev_stop_pres=stop_pres;
                
            end
            
           
        end
        
       
        
        
    end
    
    if isprint==0&icy>0
        strto_write
        fprintf(fw1,'%s\n',strto_write)
    end
    
end
disp(' ')
if ifl>0
    fprintf(fw1,'%s\n', ['\nobreakhline']);
    
end
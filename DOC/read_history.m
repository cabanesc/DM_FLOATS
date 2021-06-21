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


if Dim.n_prof.dimlength>1
 %  disp('Multi profile file')
end
ifl=0;
for n_prof=1:Dim.n_prof.dimlength
    profile_number=n_prof;
    icy=0;
    % binary ID for real time QC tests (reference table 11, Argo user manual)
    binID=[2;4;8;16;32;64;128;256;512;1024;2048;4096;8192;16384;32768;65536;131072;261144;524288;1044576;2097152;4194304];

    % name of the real time QC tests (reference table 11, Argo user manual)
    RQC_name={'Platform identification test';'Impossible Date Test';'Impossible Location Test';'Position on Land test';'Impossible Speed test';'Global Range test';'Regional Global Parameter test';'Pressure increasing test';'Spike test';'Top and bottom spike test (obsolete)';'Gradient test';'Digit rollover test';'Stuck value test';'Density inversion test';'Grey list test'; 'Gross salinity or temperature sensor drift test';'Visual qc test';'Frozen profile test';'Deepest pressure test';'Questionable Argos position test';'Near-surface unpumped CTD salinity test'; 'Near-surface mixed air/water test'};
    %  disp(' ')
    %  disp(' ')
    %  disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    %  disp(['PLATFORM NUMBER: ' strtrim(Co.platform_number.data(profile_number,:))])
    %  disp(['CYCLE NUMBER: ' num2str(Co.cycle_number.data(profile_number)) ' ' Co.direction.data(profile_number)])
    %  disp(['PROFILE NUMBER: ' num2str(k) ' ' Co.vertical_sampling_scheme.data(profile_number,:)])
    %  disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    %  disp(' ')
    %  disp(' ')

    for k=1:Dim.n_history.dimlength
        
        % HISTORY_DATE
       % keyboard
        thedate = (squeeze(Co.history_date.data(k,profile_number,:)))';
        thedatenum(k) = datenum(thedate,'yyyymmddHHMMSS');
        
    end
    [sortdate,isort]=sort(thedatenum);

    for k=[isort]
        
        
        
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
%          % HISTORY QC_TEST
%          total = hex2dec(squeeze(Co.history_qctest.data(k,profile_number,:))');
%          totala= total;
%          % find the real time tests corresponding to the hex code
%          add=[];
%          binID2=binID;
%          %if strcmp(software,'COAR (real time tests) ')|strcmp(software,'COQC')
%              
%              ikl=1;
%              %keyboard
%              while total>0
%                  % trouve le binID le plus proche
%                  [soustotal,isort] = sort(total-binID2);
%                  ik = find(soustotal<0);
%                  soustotal(ik)=NaN;
%                  [mins,is]=min(soustotal);
%                  % index des tests concernes
%                  add(ikl)=isort(is);
%                  total=total-binID(isort(is));
%                  %binID2(isort(is))=NaN;
%                  ikl=ikl+1;
%              end
%          %end
        % HISTORY_ACTION
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
        
        
        
        % PRINT ON SCREEN
        
        
        
    %      disp('---------------------------------------------------')
    %      disp(['On ' thedatestr ': ' action  ' using the software: ' software] )
    %      
        cycle = [sprintf('%3.3i',Co.cycle_number.data(profile_number,:)), strtrim(Co.direction.data(profile_number,:))];
        if isempty(previous_value)==0&strcmp(strtrim(software),'SCOOP2')
          keyboard
           % recherche le nouveau qc
           iniv=find(Co.pres.data(profile_number,:)==start_pres);
           new_flag=Co.([strtrim(lower(parameter)) '_qc']).data(profile_number,iniv);
            if str2num(new_flag)~=3
            if ifl==0
            fprintf(fw1,'%s\n', [strtrim(Co.platform_number.data(profile_number,:)) ' & ' cycle '&' strtrim(parameter) '&' num2str(previous_value) '&' new_flag '&' num2str(start_pres) ' : ' num2str(stop_pres) '&' thedatestr '\\']);
            else
              if icy==0
                 fprintf(fw1,'%s\n', [ ' & ' cycle '&' strtrim(parameter) '&' num2str(previous_value) '&' new_flag '&' num2str(start_pres) ' : ' num2str(stop_pres)  '&' thedatestr '\\']);
              else
                 fprintf(fw1,'%s\n', [ ' & &' strtrim(parameter) '&' num2str(previous_value) '&' new_flag '&' num2str(start_pres) ' : ' num2str(stop_pres)  '&' thedatestr '\\']);
              end
            end
            icy=icy+1;
            ifl=ifl+1;
            end
        % disp([' Action on: ' strtrim(parameter) ', Old flag value: ' num2str(previous_value) ', at pressures: ' num2str(start_pres) ':' num2str(stop_pres) ',' thedatestr ',' software] )
        end
        
%          %disp(totala)
%          if strcmp(action, ' THESE TESTS FAILED ');
%          for l=1:length(add)
%          if strcmp(RQC_name{add(l)},'Near-surface mixed air/water test')==0
%              disp(['The ' RQC_name{add(l)} '  failed cycle ' cycle])
%              end
%          end
%          end
        
    %    disp(' ')
    
        
    end

end
disp(' ')
if ifl>0
fprintf(fw1,'%s\n', ['\hline']);

end
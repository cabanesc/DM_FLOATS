%--------------------------------------------------------------------------
% Auteur: DAVID Nicolas
% Date: 02/05/2007
% Modification: 
% Objectif: Interpoller les données Ovide de la même manière que les
% données Argo
%--------------------------------------------------------------------------
function [CTDi] = interpolation_m(CTD,F,y_name,PARAM)	

newy=F.(y_name).data; % nbcycles x nlevels

field_names=fieldnames(CTD);
% output initialization 
for j=1:length(field_names)
	CTDi.(field_names{j})=F.(field_names{j});
	if size(CTD.(field_names{j}).data,2)==size(CTD.(y_name).data,2)
		for k=1:size(F.(y_name).data,1) % loop on nbcycles
		   CTDi.(field_names{j}).data(k,:)=[NaN]	;
		end
	else
	for k=1:size(F.(y_name).data,1) % loop on nbcycles
	CTDi.(field_names{j}).data(k,:)=CTD.(field_names{j}).data;
	end
	end
end

if strcmp(y_name,'pres')
	for j=1:length(field_names)
		for k=1:size(F.(y_name).data,1) % loop on nbcycles
		
			pression_argo = F.(y_name).data(k,:);
			pression_ctd = CTD.(y_name).data;
			[nz_ctd] = size(pression_ctd,2);
			[ncy,nz] = size(pression_argo);
			test_NaN = isfinite(pression_argo);
		
		
			if size(CTD.(field_names{j}).data,2)==nz_ctd % interpolation
				%CTDi.(field_names{j}).data(k,:)=[NaN]	;
				for i = 2:(nz-1)
					if (test_NaN(i) == 0)
						val= NaN;
					else
						clear A;
						A = [];
						[C,ind] = min(abs(pression_ctd - pression_argo(i)));
						
						if (pression_ctd(ind) < PARAM.STEP_SURF)
							A = [A,NaN];
						elseif(pression_argo(i) < PARAM.LIM_SURF_DEEP)
						
							for ll = max((ind - floor(PARAM.STEP_SURF/2)+1),1):min((ind + floor(PARAM.STEP_SURF/2)+1),nz_ctd)  
								A = [A,CTD.(field_names{j}).data(1,ll)];
							end
							val = mean(A);
						else
							for ll = max((ind - floor(PARAM.STEP_DEEP/2)+1),1):min((ind + floor(PARAM.STEP_DEEP/2)+1),nz_ctd) 
								A = [A,CTD.(field_names{j}).data(1,ll)];
							end
							val = mean(A);
						end
						val = mean(A);
					end
					
					CTDi.(field_names{j}).data(k,i) = val;
				end
				CTDi.(field_names{j}).data(k,i) = NaN;
			% else
				% CTDi.(field_names{j}).data(k,:)=CTD.(field_names{j}).data;
			end
			CTDi.(y_name).data(k,:)=F.(y_name).data(k,:);
				
			
		end
	end
elseif strcmp(y_name,'ptmp0')|strcmp(y_name,'ptmp1000')
disp('interp_climatology')
	for k=1:size(F.(y_name).data,1) % loop on nbcycles	
        CTDi.pres.data(k,:)=[NaN];    
		CTDi.psal.data(k,:)=[NaN]; 
		[S_h,P_h]=interp_climatology(CTD.psal.data(1,:)',CTD.(y_name).data(1,:)',CTD.pres.data(1,:)',F.psal.data(k,:)',F.(y_name).data(k,:)',F.pres.data(k,:)'); % 
		CTDi.psal.data(k,:)=S_h';
		CTDi.pres.data(k,:)=P_h';
		CTDi.(y_name).data(k,:)=F.(y_name).data(k,:);
	end
	   
else
	error (['interpolation no implemented for ' y_name])
end





function [ pt_input_parameters ] = load_configuration( ps_filename , flt_name, DIR_DATA, VERSION_OW, DIR_OWC)

%
% loadConfiguration
%
% loadConfiguration loads a set of properties from a file
% and returns the results in a structure variable.
%
% The file is formatted as value pairs (one line per value):
%
%    key_name=key_value
%
% Also, any line starting with a '%' will be ignored.
%
% Note: all values of key_name must be valid matlab
%       variable names.
%
%       In the case of duplicate entries, the last key
%       to appear in the file will override previous
%       definitions.
%
%       Invalid lines are disregarded.
%
%  Jason Fabritz, 2001
%
% 16/10/2007 V. Thierry 
% ajout d'un argument d'entree 'flt_name' et chnagement du nom du reportoir des figures
% il dépend maintenant du numéro WMO du flotteur.
% 15/10/2012 C.Cabanes
% ajout de 2 arguments d'entree DIR_DATA (repertoire des Donnees/Resultats) et VERSION_OW (indique la version OW utilisée
% 21.09.2020 ==> DIR_OWC added by T. Reynaud

if nargin==2
   DIR_DATA='';
   VERSION_OW='';
end
if nargin==3
   VERSION_OW='';
end
if nargin==4
   DIR_OWC='';
end


lh_file = fopen( ps_filename ) ;

pt_input_parameters = struct( 'CONFIGURATION_FILE', ps_filename ) ;

while not( feof( lh_file ) ) ;   
   ls_line = fgetl( lh_file ) ;
   if length( ls_line ) > 0 ;
      if ls_line( 1:1 ) ~= '%' ;
	   	ln_equals_index = findstr( '=', ls_line ) ;
		   if length( ln_equals_index ) > 0  ;
      		if ln_equals_index( 1 ) > 1 ;
   	   	        sKey = ls_line( 1 : ln_equals_index( 1 ) - 1 ) ;
   	   	        
   	   	
   	   	   
                %VT
                if strcmp(sKey,'FLOAT_PLOTS_DIRECTORY')
                    sValue = [deblank( ls_line( ln_equals_index( 1 ) + 1 : length( ls_line ) ) ) flt_name '/'] ;
                else
	         	  sValue = deblank( ls_line( ln_equals_index( 1 ) + 1 : length( ls_line ) ) ) ;
                end
                %VT
                
                %CC
   	   	if isempty(strfind(sValue,'{DIR_DATA}'))==0
   	   	    sValue=strrep(sValue,'{DIR_DATA}',DIR_DATA);
   	   	end
   	   	
   	   	if isempty(strfind(sValue,'{VERSION_OW}'))==0
   	   	    sValue=strrep(sValue,'{VERSION_OW}',VERSION_OW);
   	   	end
   	   	%CC
        
        % Added by T. Reynaud 21.09.2020
   	   	if isempty(strfind(sValue,'{DIR_OWC}'))==0
   	   	    sValue=strrep(sValue,'{DIR_OWC}',DIR_OWC);
   	   	end
                
		         pt_input_parameters = setfield( pt_input_parameters, sKey, sValue ) ;
	      	end
   		end
   	end
   end   
end


fclose( lh_file ) ;

return


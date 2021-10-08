% genere Configration table tex
function generate_conf_table(numconfig,refdatabase)
%init_path
C = load_configuration('config.txt');
DIR_FTP=C.DIR_FTP;
DIR_PLOT=C.DIR_PLOT;
DIR_DATA=C.DIR_DATA;
DIR_PLOTOW=[C.DIR_DATA 'float_plots/CONFIG'];
DIR_DMQC=C.DIR_DM_FILES;
DIR_CODES=C.DIR_CODES;
lo_system_configuration = load_configuration( [DIR_CODES 'LPO_CODES_ATLN_NEW/ow_config/ow_config_' num2str(numconfig) '.txt'] );

fw1=fopen(['./REPORTS/table_config' num2str(numconfig) '.tex'],'w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TABLE3
fprintf(fw1,'%s\n', ['\begin{table}[h]']);
fprintf(fw1,'%s\n', ['$$']);
fprintf(fw1,'%s\n', ['\begin{tabular}{|l|c|}']);
fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['OW CONFIGURATION 		& ' num2str(numconfig) '     	\\']);
fprintf(fw1,'%s\n', ['				&                \\']);
fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['CONFIG\_MAX\_CASTS		& ' lo_system_configuration.CONFIG_MAX_CASTS '     	\\']);
fprintf(fw1,'%s\n', ['MAP\_USE\_PV			& ' lo_system_configuration.MAP_USE_PV '       	\\']);
fprintf(fw1,'%s\n', ['MAP\_USE\_SAF		        & ' lo_system_configuration.MAP_USE_SAF '        	\\']);
fprintf(fw1,'%s\n', ['MAPSCALE\_LONGITUDE\_LARGE	& ' lo_system_configuration.MAPSCALE_LONGITUDE_LARGE '     	\\']);
fprintf(fw1,'%s\n', ['MAPSCALE\_LONGITUDE\_SMALL	& ' lo_system_configuration.MAPSCALE_LONGITUDE_SMALL '        \\']);
fprintf(fw1,'%s\n', ['MAPSCALE\_LATITUDE\_LARGE 	& ' lo_system_configuration.MAPSCALE_LATITUDE_LARGE '           \\']);
fprintf(fw1,'%s\n', ['MAPSCALE\_LATITUDE\_SMALL 	& ' lo_system_configuration.MAPSCALE_LATITUDE_SMALL '      \\']);
if lo_system_configuration.MAP_USE_PV==1
fprintf(fw1,'%s\n', ['MAPSCALE\_PHI\_LARGE	 	& ' lo_system_configuration.MAPSCALE_PHI_LARGE '      \\']);
fprintf(fw1,'%s\n', ['MAPSCALE\_PHI\_SMALL	 	& ' lo_system_configuration.MAPSCALE_PHI_SMALL '    \\']);
end
fprintf(fw1,'%s\n', ['MAPSCALE\_AGE		 	& ' lo_system_configuration.MAPSCALE_AGE '    \\']);
fprintf(fw1,'%s\n', ['MAPSCALE\_AGE\_LARGE		& ' lo_system_configuration.MAPSCALE_AGE_LARGE '    	\\']);
fprintf(fw1,'%s\n', ['MAP\_P\_EXCLUDE		 	& ' lo_system_configuration.MAP_P_EXCLUDE '      \\']);
fprintf(fw1,'%s\n', ['MAP\_P\_DELTA		 	& ' lo_system_configuration.MAP_P_DELTA '      \\']);
if isempty(strfind(refdatabase,'&'))
fprintf(fw1,'%s\n', ['Reference data base      	&  ' refdatabase '   \\ ']);
else
    [a,b]=strtok(refdatabase,'\&');
 fprintf(fw1,'%s\n', ['Reference data base      	&  ' a '   \\ ']);  
 fprintf(fw1,'%s\n', ['                         	&  ' b '   \\ ']); 
end
 
fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['\end{tabular}']);
fprintf(fw1,'%s\n', ['$$']);

%fprintf(fw1,'%s\n', ['\caption{Parameters of the OWC method }']);
fprintf(fw1,'%s\n', ['\label{tab3}']);
fprintf(fw1,'%s\n', ['\end{table}']);

fclose(fw1);

fw2=fopen(['./REPORTS/ref_config' num2str(numconfig) '.txt'],'w');
str=strrep(refdatabase,'\&','&');
fprintf(fw1,'%s\n',str);
fclose(fw2)

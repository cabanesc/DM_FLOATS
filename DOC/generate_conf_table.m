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

fprintf(fw1,'%s\n', ['\documentclass[11pt,titlepage]{article}']);
fprintf(fw1,'%s\n', ['\usepackage[latin1]{inputenc}']);
%\usepackage[french]{babel}
fprintf(fw1,'%s\n', ['\usepackage[dvips,final]{graphicx}']);

fprintf(fw1,'%s\n', ['\usepackage{color}']);
fprintf(fw1,'%s\n', ['\usepackage{tabularx}']);
fprintf(fw1,'%s\n', ['\usepackage{float}']);
%\usepackage{rotfloat}
fprintf(fw1,'%s\n', ['\usepackage{calc}']);
fprintf(fw1,'%s\n', ['\usepackage{soul}']);
fprintf(fw1,'%s\n', ['\usepackage[small,hang]{caption}']);
fprintf(fw1,'%s\n', ['\usepackage{subcaption}']);
fprintf(fw1,'%s\n', ['\usepackage[colorlinks={true}]{hyperref}']);
fprintf(fw1,'%s\n', [' ']);
fprintf(fw1,'%s\n', ['\setlength{\textwidth}{430pt}']);
fprintf(fw1,'%s\n', ['\setlength{\oddsidemargin}{0pt}']);
fprintf(fw1,'%s\n', ['\setlength{\evensidemargin}{0pt}']);
fprintf(fw1,'%s\n', ['\setlength{\marginparwidth}{0pt}']);
%\setlength{\footskip}{20pt}
fprintf(fw1,'%s\n', ['\setlength{\headsep}{35pt}']);
fprintf(fw1,'%s\n', ['\setlength{\topmargin}{-15pt}']);
fprintf(fw1,'%s\n', ['\setlength{\textheight}{650pt}']);

fprintf(fw1,'%s\n', ['\newcommand{\dg}{$^\circ$}']);
fprintf(fw1,'%s\n', ['\newcommand{\Sv}{$\times 10^6$ m$^3$.s$^{-1}$}']);
fprintf(fw1,'%s\n', ['\newcommand{\pw}{$\times 10^{15}$ W}']);
fprintf(fw1,'%s\n', ['\newcommand{\cmps}{cm.s$^{-1}$}']);
fprintf(fw1,'%s\n', ['\newcommand{\mps}{m.s$^{-1}$}']);
fprintf(fw1,'%s\n', ['\newcommand{\mcps}{m$^2$.s$^{-1}$}']);
fprintf(fw1,'%s\n', ['\newcommand{\upv}{$\times 10^{-11}$ m$^{-1}$.s$^{-1}$}']);
%\newcommand{\sig0}{$\sigma_0$}
fprintf(fw1,'%s\n', ['\newcommand{\usig}{kg.m$^{-3}$}']);


%\graphicspath{{/home1/corsen/perso/ccabanes/dvlpRD/Argo/TD/OW/data/float_plots/CONFIG1/1900075/}  {/home1/corsen/perso/ccabanes/Results/Argo/TD/Traitements/Plot/Verif_flag/1900075/}}

 fprintf(fw1,'%s\n', ['\usecounter{page}']);
 fprintf(fw1,'%s\n', ['\begin{document}']);

fprintf(fw1,'%s\n', ['\pagestyle{plain}']);
fprintf(fw1,'%s\n', ['\pagenumbering{arabic}']);

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

fprintf(fw1,'%s\n', ['\end{document}']);

fclose(fw1);

fw2=fopen(['./REPORTS/ref_config' num2str(numconfig) '.txt'],'w');
str=strrep(refdatabase,'\&','&');
fprintf(fw1,'%s\n',str);
fclose(fw2)

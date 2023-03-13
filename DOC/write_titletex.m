function write_titletex(fw1,TITRE,DOI_ARGO,AUTHORS,file_cor,float_list)

fprintf(fw1,'%s\n', ['\documentclass[11pt,titlepage]{article}']);
fprintf(fw1,'%s\n', ['\usepackage[latin1]{inputenc}']);
%\usepackage[french]{babel}
fprintf(fw1,'%s\n', ['\usepackage[dvips,final]{graphicx}']);

fprintf(fw1,'%s\n', ['\usepackage{color}']);
fprintf(fw1,'%s\n', ['\usepackage{tabularx}']);
fprintf(fw1,'%s\n', ['\usepackage{longtable}'])
fprintf(fw1,'%s\n', ['\usepackage{float}']);
%\usepackage{rotfloat}
fprintf(fw1,'%s\n', ['\usepackage{calc}']);
fprintf(fw1,'%s\n', ['\usepackage{soul}']);
fprintf(fw1,'%s\n', ['\usepackage[small,hang]{caption}']);
fprintf(fw1,'%s\n', ['\usepackage{subcaption}']);
fprintf(fw1,'%s\n', ['\usepackage[colorlinks={true}]{hyperref}']);

fprintf(fw1,'%s\n', [' ']);
fprintf(fw1,'%s\n', ['\makeatletter']);
fprintf(fw1,'%s\n', ['\def\nobreakhline{\multispan\LT@cols\unskip\leaders\hrule\@height\arrayrulewidth\hfill\cr\noalign{\penalty10000}}']);
fprintf(fw1,'%s\n', ['\makeatother']);
fprintf(fw1,'%s\n', [' ']);

fprintf(fw1,'%s\n', ['\setlength{\textwidth}{430pt}']);
fprintf(fw1,'%s\n', ['\setlength{\oddsidemargin}{0pt}']);
fprintf(fw1,'%s\n', ['\setlength{\evensidemargin}{0pt}']);
fprintf(fw1,'%s\n', ['\setlength{\marginparwidth}{0pt}']);
%\setlength{\footskip}{20pt}
fprintf(fw1,'%s\n', ['\setlength{\headsep}{35pt}']);
fprintf(fw1,'%s\n', ['\setlength{\topmargin}{-15pt}']);
fprintf(fw1,'%s\n', ['\setlength{\textheight}{650pt}']);
%\setlength{\headheight}{10pt}

% Pour la page de garde
%\newlength{\header}
%\setlength{\header}{\headheight + \headsep}

%\renewcommand{\captionfont}{\footnotesize}
%\renewcommand{\floatpagefraction}{0.9}
%\renewcommand{\textfraction}{.1}

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
%\setcounter{page}{37}
% d\'ebut page de garde
fprintf(fw1,'%s\n', ['\begin{titlepage}']);
%\thispagestyle{empty}
%\vspace*{-\header}

fprintf(fw1,'%s\n', ['\begin{center}']);
fprintf(fw1,'%s\n', ['{']);
fprintf(fw1,'%s\n', ['\newfont{\gtitre}{cmssbx10 at 16pt}']);
fprintf(fw1,'%s\n', ['\newfont{\mtitre}{cmss10 at 14pt}']);

fprintf(fw1,'%s\n', ['\begin{tabularx}{\textwidth}{lXr}']);

fprintf(fw1,'%s\n', ['\vspace*{2cm} & & \\']);
fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['\vspace*{0.5cm} & & \\']);
fprintf(fw1,'%s\n', ['\multicolumn{3}{c}{\gtitre ' TITRE  ' } \\']);
fprintf(fw1,'%s\n', ['\vspace*{-0.2cm} & & \\']);
fprintf(fw1,'%s\n', ['\multicolumn{3}{c}{\gtitre  ' DOI_ARGO '} \\']);
fprintf(fw1,'%s\n', ['\vspace*{0.5cm} & & \\']);
fprintf(fw1,'%s\n', ['\hline']);
fprintf(fw1,'%s\n', ['\vspace*{0.5cm} & & \\']);
fprintf(fw1,'%s\n', ['\multicolumn{3}{c}{ ' AUTHORS ' }\\ ']);
fprintf(fw1,'%s\n', ['\vspace*{0.5cm} & & \\']);
fprintf(fw1,'%s\n', ['\multicolumn{3}{c}{ SO ARGO - LOPS report - Update \today}\\  ']);
fprintf(fw1,'%s\n', ['\end{tabularx}']);
fprintf(fw1,'%s\n', [' ']);
fprintf(fw1,'%s\n', ['\vspace*{2cm} ']);
fprintf(fw1,'%s\n', [' \textbf{Summary}']);


if exist(file_cor)~=0
    [num,wmo_corr,col1,col2_1,col2_2,col3,col4,col5]=get_txtfile_col(file_cor,';');

    if sum(ismember(wmo_corr,float_list))>=1 % on fait un tableau
        % preparation du tableau
        fprintf(fw1,'%s\n', ['\renewcommand\arraystretch{1.5}']);
        fprintf(fw1,'%s\n', ['\begin{table}[h]']);
        fprintf(fw1,'%s\n', ['$$']);
        fprintf(fw1,'%s\n', ['\begin{tabular}{|l|c|}']);
        fprintf(fw1,'%s\n', ['\hline']);
        fprintf(fw1,'%s\n', ['WMO Number & DM Salinity Correction \\']);
        fprintf(fw1,'%s\n', ['\hline']);
        fprintf(fw1,'%s\n', ['\hline']);
        %--------------------
        % remplissage du tableau
        for ik=1:length(float_list)
            iil=find(findstr_tab(wmo_corr,float_list{ik}));
            if isempty(iil)==0
                fprintf(fw1,'%s\n',[float_list{ik} ' & ' col4{iil} ' \\']);
                fprintf(fw1,'%s\n', ['\hline']);
            end
        end
        % fin du tableau et legende
        fprintf(fw1,'%s\n', ['\end{tabular}']);
        fprintf(fw1,'%s\n', ['$$']);
        fprintf(fw1,'%s\n', ['\caption{Salinity Correction applied in delayed mode for each float.}']);
        fprintf(fw1,'%s\n', ['\label{tab0}']);
        fprintf(fw1,'%s\n', ['\end{table}']);
        fprintf(fw1,'%s\n', ['\renewcommand\arraystretch{1.2}']);
    end
end

% \begin{figure}[h]
% \begin{center}
% $$
% \includegraphics[width=10cm]{1900383_pos_flagused.eps}
% $$
% \end{center}
% \end{figure}

fprintf(fw1,'%s\n', ['}']);
fprintf(fw1,'%s\n', ['\end{center}']);


fprintf(fw1,'%s\n', ['\end{titlepage}']);

fprintf(fw1,'%s\n', ['\tableofcontents']);
fprintf(fw1,'%s\n', ['\clearpage']);
fprintf(fw1,'%s\n', ['\setcounter{section}{0}']);
function init_path(action,rep,rdir)

%path(pathdef)

switch lower(rep)
    case {lower('CHANGE_FLAGS')}
        switch lower(action)
            case lower('add')
                close all
                addpath(fullfile(rdir,'/lib/seawater_330_its90'));
                addpath(fullfile(rdir,'/lib/libargo'));
                addpath(fullfile(rdir,'/lib/ext_lib'));
                addpath(fullfile(rdir,'CHANGE_FLAGS'));
                addpath(fullfile(rdir,'CHANGE_FLAGS/util'));
                addpath(fullfile(rdir));
                addpath('/home1/homedir5/perso/ccabanes/matlab/util/m_map1.4m');
            case lower('clear')
                path(pathdef)
        end
        
    case {lower('PROG_QC2015')}
        switch lower(action)
            case lower('add')
                addpath(fullfile(rdir,'/lib/seawater_330_its90'));
                addpath(fullfile(rdir,'/lib/libargo'));
                addpath(fullfile(rdir,'/lib/ext_lib'));
                addpath(fullfile(rdir,'/PROG_QC2015/util'));
                addpath(fullfile(rdir,'/lib/tool_Ocean/cmocean'));
                addpath('/home1/homedir5/perso/ccabanes/matlab/util/m_map1.4m');
            case lower('clear')
                path(pathdef)
        end
    case {lower('VERIF_FLAG')}
        switch lower(action)
            case lower('add')
                addpath(fullfile(rdir,'/lib/seawater_330_its90'));
                addpath(fullfile(rdir,'/lib/libargo'));
                addpath(fullfile(rdir,'/lib/ext_lib'));
                addpath(fullfile(rdir));
                addpath(fullfile(rdir,'/VERIF_FLAG/util'));
                addpath(fullfile(rdir,'/lib/tool_Ocean/cmocean'));
                addpath('/home1/homedir5/perso/ccabanes/matlab/lib_lpo/hydcal');
                addpath('/home1/homedir5/perso/ccabanes/matlab/util/m_map1.4m');
            case lower('clear')
                path(pathdef)
        end
    case {lower('VERIF_PROF1')}
        switch lower(action)
            case lower('add')
                addpath(fullfile(rdir,'/lib/seawater_330_its90'));
                addpath(fullfile(rdir,'/lib/libargo'));
                addpath(fullfile(rdir,'/lib/ext_lib'));
                addpath(fullfile(rdir));
                addpath(fullfile(rdir,'/VERIF_PROF1/util'));
                addpath(fullfile(rdir,'/VERIF_PROF1/config_campagne'));
                addpath(fullfile(rdir,'/lib/tool_Ocean/cmocean'));
                addpath(fullfile(rdir,'/lib/gsw_matlab_v3_01'));
                addpath('/home1/homedir5/perso/ccabanes/matlab/util/m_map1.4m');
            case lower('clear')
                path(pathdef)
        end
    case {lower('CPCOR')}
        switch lower(action)
            case lower('add')
                addpath(fullfile(rdir,'/lib/seawater_330_its90'));
                addpath(fullfile(rdir,'/lib/libargo'));
                addpath(fullfile(rdir,'/lib/ext_lib'));
                addpath(fullfile(rdir));
                addpath(fullfile(rdir,'/CPCOR/util'));
                addpath(fullfile(rdir,'/VERIF_PROF1/config_campagne'));
                addpath(fullfile(rdir,'/lib/tool_Ocean/cmocean'));
				addpath(fullfile(rdir,'/lib/gsw_matlab_v3_01'));
                addpath('/home1/homedir5/perso/ccabanes/matlab/util/m_map1.4m');
            case lower('clear')
                path(pathdef)
        end
      case {lower('COMPARE_FLOAT_REF_TR')}
        switch lower(action)
            case lower('add')
                addpath(fullfile(rdir,'/lib/seawater_330_its90'));
                addpath(fullfile(rdir,'/lib/libargo'));
                addpath(fullfile(rdir,'/lib/ext_lib'));
                addpath(fullfile(rdir));
                addpath(fullfile(rdir,'/COMPARE_FLOAT_REF_TR/util'));
            case lower('clear')
                path(pathdef)
        end  
     case {lower('LPO_CODES_ATLN_NEW')}
        switch lower(action)
            case lower('add')
                close all
                addpath(fullfile(rdir,'/lib/seawater_330_its90'));
                addpath(fullfile(rdir,'/lib/libargo'));
                addpath(fullfile(rdir,'/lib/ext_lib'));
                addpath(fullfile(rdir,'LPO_CODES_ATLN_NEW'));
                addpath(fullfile(rdir,'LPO_CODES_ATLN_NEW/util'));
                addpath(fullfile(rdir));
                addpath('/home1/homedir5/perso/ccabanes/matlab/util/m_map1.4m');
            case lower('clear')
                path(pathdef)
        end
    case {lower('DOC')}
        switch lower(action)
            case lower('add')
                addpath(fullfile(rdir,'/lib/seawater_330_its90'));
                addpath(fullfile(rdir,'/lib/libargo'));
                addpath(fullfile(rdir,'/lib/ext_lib'));
                addpath(fullfile(rdir));
                addpath(fullfile(rdir,'/lib/tool_Ocean/cmocean'));
                addpath('/home1/homedir5/perso/ccabanes/matlab/util/m_map1.4m');
                addpath(fullfile(rdir,'/PROG_QC2015'));
            case lower('clear')
                path(pathdef)
        end        
    case {lower('CORRECTIONS')}
        switch lower(action)
            case lower('add')
                addpath(fullfile(rdir,'/CORRECTIONS/lib/seawater_330_its90'));
                addpath(fullfile(rdir,'/CORRECTIONS/lib/'));
                addpath(fullfile(rdir));
            case lower('clear')
                path(pathdef)
        end
end

eval(['cd ',rdir]);

return



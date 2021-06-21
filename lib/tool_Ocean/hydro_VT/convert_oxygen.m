% fonction qui convertit les données d'oxygène d'une unité à l'autre
% function [valout]=convert_oxygen(valin,uin,uout,varargin)
% valout = valeur d'oxygene converties dans l'unité uout
% valin = valeur d'oxygene en entrée dans l'unité uin
% varargin: anomalie de densité potentielle référencée à 0 de l'eau considérée, si aucun valeur n'est entrée, on
% considère qu'on a de l'eau douce et l'anomalie de densité est nulle
% (rho=1000) sinon rho =1000+varargin
% uin et uout sont:
% mL/L: milli litre par litre
% mmol/m3: milli mole par mètre cube
% mumol/L: micromole par litre
% mg/L : milli gramme par litre
% mumol/kg : micromole par kilo
% mL/L * 44.66 = mmol/m3 = mumol/L
% mL/L*1.42903 = mg/L
% mg/L * 44.66/1.42903 = mumol/L =mmol/m3
% mumol/L = (rho/1000)* mumol/kg 



function [valout]=convert_oxygen(valin,uin,uout,varargin)

coef1=44.66;
coef2=1.42903;
if nargin == 3
    rho=1000;
elseif nargin == 4
    rho=varargin{1}+1000;
end
coef3=1000/rho;

valout=[];

switch uin
    case 'mL/L'
        switch uout
            case {'mmol/m3','mumol/L'}
                valout=valin*coef1;
            case 'mg/L'
                valout=valin*coef2;
            case 'mumol/kg'
                valout=valin*coef1*coef3;
        end
    case {'mmol/m3','mumol/L'}
        switch uout
            case 'mL/L'
                valout=valin/coef1;
            case 'mg/L'
                valout=valin*coef2/coef1;
            case 'mumol/kg'
                valout=valin*coef3;
        end
    case 'mg/L'
        switch uout
            case 'mL/L'
                valout=valin/coef2;
            case {'mmol/m3','mumol/L'}
                valout=valin*coef1/coef2;
            case 'mumol/kg'
                valout=valin*coef3*coef1/coef2;
        end
    case 'mumol/kg'
        switch uout
            case 'mL/L'
                valout=valin/(coef1*coef3);
            case {'mmol/m3','mumol/L'}
                valout=valin/coef3;
            case 'mg/L'
                valout=valin*coef2/(coef3*coef1);
        end
end
           
 



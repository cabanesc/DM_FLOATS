% fonction qui convertit les donn�es d'oxyg�ne d'une unit� � l'autre
% function [valout]=convert_oxygen(valin,uin,uout,varargin)
% valout = valeur d'oxygene converties dans l'unit� uout
% valin = valeur d'oxygene en entr�e dans l'unit� uin
% varargin: anomalie de densit� potentielle r�f�renc�e � 0 de l'eau consid�r�e, si aucun valeur n'est entr�e, on
% consid�re qu'on a de l'eau douce et l'anomalie de densit� est nulle
% (rho=1000) sinon rho =1000+varargin
% uin et uout sont:
% mL/L: milli litre par litre
% mmol/m3: milli mole par m�tre cube
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
           
 



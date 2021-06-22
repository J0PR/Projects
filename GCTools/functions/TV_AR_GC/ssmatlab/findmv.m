function [J,miss]=findmv(y)
%
% this function finds missing values in vector y.
% Missing values in y are entered as NaN.
%
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos, 
% Subdireccion Gral. de Analisis y P.E., 
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should 
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
miss=0;
p=length(y);
JJ=eye(p);
J=[];
for i=1:p
 if isnan(y(i))
  miss=miss+1;
 else
  J=[J;JJ(i,:)];
 end
end
y(isnan(y))=0;

function FTsh=shiftFT(FT,vect) ;
%  FTsh=shiftFT(FT,vect) ;
% desplaza todos los coeficientes de la transformada de Fourier.
% pone el coeficiente (1,1) en la posicion vect(1:2)+(1,1)
% similar a la funccion fftshift
%
% FT : fourier transform a desplazar
% vect = vector de desplazamiento
% FTsh = resultado
%
%Sylvain Fischer,  2000


tam= size(FT) ;
vect=mod(vect,tam)  ;
FTsh=FT ;

FTsh(1:vect(1),1:vect(2))=FT(tam(1)-vect(1)+1:tam(1),tam(2)-vect(2)+1:tam(2)) ;

FTsh(1+vect(1):tam(1),1:vect(2))=FT(1:tam(1)-vect(1),tam(2)-vect(2)+1:tam(2)) ;

FTsh(1:vect(1),1+vect(2):tam(2))=FT(tam(1)-vect(1)+1:tam(1),1:tam(2)-vect(2)) ;

FTsh(1+vect(1):tam(1),1+vect(2):tam(2))=FT(1:tam(1)-vect(1),1:tam(2)-vect(2)) ;

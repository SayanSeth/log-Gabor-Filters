function invg=symetriefft(g,opt) ;
%invg=symetriefft(g,opt) ;
% devolve el symetrico de la imagen fourier g
%el zero del espacio esta cerca del centro de la imagen
%
%opt=0 : invg(1,:)=invg(:,1)=0
%opt=1 : invg(1,:)=conj(g(1,:)) ...
%
%sylvain 16-10-2003 / 11-04-2004
%funccion a mejorar

tam=size(g) ;
invg=zeros(tam) ;
if exist('opt')==1,
else opt=0 ;
end

if mod(tam(1),2)==0, %imagen par
    invg(2:tam(1),2:tam(2))=conj(g(tam(1):-1:2,tam(2):-1:2)) ;     
    if opt==1,
         invg(1,2:tam(2))=conj(g(1,tam(2):-1:2)) ;
         invg(2:tam(1),1)=conj(g(tam(1):-1:2,1)) ; 
         invg(1,1)=conj(g(1,1)) ;
    end
else %impar
    invg(1:tam(1),1:tam(2))=conj(g(tam(1):-1:1,tam(2):-1:1)) ;     
end


function [coefG,csf] = dcmp(img,csf,coefRose,divD,roseta)

% function [coefG,csf] = dcmp(img,csf,coefRose,divD,roseta)
%
%   This function DECOMPOSES an image into a log-Gabor subbands.
%
% INPUT
%       img - input image
%       csf - contrast sensitivity function (-1 for none)
%       coefRose - log-Gabor coefficients arranged as logGaborFilters.m
%       divD - normalization factor (recommended =0)
%       roseta - log-Gabor filter summatory delivered by logGaborFilters.m
%
% OUTPUT 
%       coefG - log-Gabor coefficients of the transform image
%       csf (additional)
%
% Note rcmp.m requires shiftFT.m
%
% Sylvain Fischer - 2007

%[coefG,csf] = dcmp(img,csf,coefRose,divD,roseta)
% procesa la transformation lineal 
% img es la imagen (en gris nxn pixels) a comprimir 
% csf : Contrast sensitivity function mxm pixels  (m=n+1 is n par, m=n si n impar)
% csf=-1 : no multiplicar por una csf
% coefRose son los filtros de la transformada wavelet produced por la function  Rosette
% mejor divD=0 (=param(17)) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                1:   parametros y opciones                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

tam= size(coefRose) ;
tam(1:2)=[real(coefRose(tam(3)+1,2,tam(3),3)) imag(coefRose(tam(3)+1,2,tam(3),3))] ;
sizefft=coefRose(tam(3)+1,4:5,tam(3),3)  ;% tamano de la imagen+border



% coordenadas de las ventanas de cada canal
v=coefRose(1:tam(3),2:2*tam(4)+1,tam(3),3)  ; 
w(:,:,1)=real(v(:,1:2:2*tam(4)-1)); w(:,:,2)=real(v(:,2:2:2*tam(4)));%Xmin max
w(:,:,3)=imag(v(:,1:2:2*tam(4)-1)); w(:,:,4)=imag(v(:,2:2:2*tam(4)));%Ymin max
W(:,:,1)= w(:,:,2)-w(:,:,1)+1 ; W(:,:,2)= w(:,:,4)-w(:,:,3) +1 ;
ampl=coefRose(tam(3)+2:2*tam(3)+1,1:tam(4),tam(3),3) ;% determinar la amplitud con la imagen gabors
if ampl==0, %no hay nada definido
    ampl(:,:)=1;
   % ampl=exp([1:tam(3)]'*log(4))*ones(1,tam(4)) ;
   % ampl(tam(3),2)=1 ;
end

coefG = coefRose ;
fftimg= fftshift(fft2(img)) ; 


if nargin>=5,
  if divD==1,   %dividir por la roseta.
    roseta(roseta==0)=0.001 ;
    fftimg =fftimg./roseta ;
  end
end

 % multiplicar por la csf
if csf~=-1 ,
  %%%%%%%%%%%%
  % paramQ global function. better to be removed before ending
  global paramQ ;
  paramQ(4)=1 ;%  umbralcsf=1 ;
  paramQ(11)=1;%  umbralJPEG=1 ;
  paramQ(12)=75;%  qualityJPEG=75 ;
  paramQ(13)=2.08;%  tasa de  compression JPEG2000=2bpp ;
  %  precision=5 ;
  paramQ(2)=32;%  global resolucion ;  resolucion=32 ;           % resolucion pantalla ,  (pixel por cm)
  paramQ(3)=40;%  global distancia;      distancia=40;         % distancia de observacion (cm)
  paramQ(1)=3;%  global correctionCSF ;   correctionCSF=3 ; %CSF con low-pass correction % 0 : sin CSF, 1 : CSF de Rust ; 2 correction paso bajo de la CSF de Rust.
  paramQ(5)=.8;%   global alfa ; alfa=.8 ;
  paramQ(6)=2;% modo de applicacion de la CSF
  paramQ(7)=1 ; % no eliminar las cadenas pequeñas en la Quantization
  paramQ(8)=3 ; % 1:quantization real/imag 2,3,4: quantiz. modulus/4,8,16,32 phases
  paramQ(15)=4 ; %  eliminar los coef aislados de Quantization step<8 ;
  paramQ(16)=2 ; %  eliminar los coef en cadenas de 2 o - de Quantization step<4 ;
  paramQ(17)=1.5 ; %  eliminar los coef en cadenas de 3 o - de Quantization step<2 ;
  paramQ(18)=0 ; %  eliminar los coef en cadenas de 4 o - de Quantization step<1 ;
  paramQ(19)=0 ; %  eliminar los coef en cadenas de 5 o - de Quantization step<1 ;
  if size(csf)~=size(fftimg) , 
        csf=csfor(tam,sizefft) ; %poner una opcion para no tener CSF
     if paramQ(6)~=1,   %csf por pesos
        for k=1:tam(3),%por promedio
          for l=1:tam(4), 
            if W(k,l,1)>1,   %si  hay filtro para este [k,l]
                csf(k,l)=sum(sum(abs(csf(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4)).*coefRose(1:W(k,l,1),1:W(k,l,2),k,l))))/sum(sum(abs(coefRose(1:W(k,l,1),1:W(k,l,2),k,l)))) ;
            end
          end
        end
      csf= csf(1:tam(3),1:tam(4));
      csf(tam(3),1)=1 ;
     end
  end
  if paramQ(6)==1, %CSF by frequency product
      fftimg=fftimg.*csf ;
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 : convolucion filtro con imagen + submuestreo para cada canal + correction de amplitud %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:tam(3)
   for l=1:tam(4) %round(1+ (tam(4)-1)*(k<tam(3)))    
       if W(k,l,1)>1,   %si  hay filtro para este [k,l]
            vect(1)=ceil((tam(1)+1)/2)-w(k,l,1) ; %shift para conservar el centro
            vect(2)=ceil((tam(2)+1)/2)-w(k,l,3) ;
            prod =fftimg(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4)).*coefRose(1:W(k,l,1),1:W(k,l,2),k,l);
            coefG(1:W(k,l,1),1:W(k,l,2),k,l)=ifft2(shiftFT(prod,-vect))*ampl(k,l);% ;%*4^(-k) ; %valores de coefG en contraste
     %       if csf~=-1,
     %           if paramQ(6)==2, %CSF por pesos
     %              coefG(:,:,k,l)=coefG(:,:,k,l); % *csf(k,l) ;%ya no se multiplica aqui por la CSF pero solo en la cuantization.
     %           end
     %       end    
       end
   end
end
%coefG(:,:,tam(3),2)=coefG(:,:,tam(3),2)*4^(tam(3)) ; % cancelar el 4^(-k) del paso alto

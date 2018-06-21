function [imgR,imgRb]= rcmp(coefG,divR,csf,roseta,coefRose,canal)

% function [imgR,imgRb] = rcmp(coefG,divR,csf,roseta,coefRose,channel)
%
%   This function RECONSTRUCTS an image from the log-Gabor transformation.
%
% INPUT
%       coefG - log-Gabor coefficients of the transform image
%       divD - normalization factor (recommended =1)
%       csf - contrast sensitivity function (-1 for none)
%       roseta - log-Gabor filter summatory delivered by logGaborFilters.m
%       coefRose - log-Gabor coefficients arranged as logGaborFilters.m
%       channel = 0 reconstrucs all subbands or [scale orient] for one subband.
%
% OUTPUT 
%       imgR - reconstructed image
%       imgRb (additional) - reconstructed image with border
%
% Note rcmp.m requires shiftFT.m
%
% Sylvain Fischer - 2007

% [imgR,imgRb] = rcmp(coefG,divR,csf,roseta,coefRose,canal)
% csf =-1:  sin dividir por la csf
% canal = 0  para reconstruir todos los canales
% canal= [escala orientacion] para recontruir 1 solo canal
% canal=( [tam(3) 2]) y [tam(3) 1] para el paso alto y el paso bajo respectivemente
%canal= -1 : no quitar el borde
%imgRb: image with the border


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          (1)            parametros                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global param;
global paramQ;
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
%parametros y variables
tam=size(coefRose) ;
tam(1:2)=[real(coefRose(tam(3)+1,2,tam(3),3)) imag(coefRose(tam(3)+1,2,tam(3),3))];
sizefft=coefRose(tam(3)+1,4:5,tam(3),3) ; % tamano de la imagen+1+border
border = coefRose(tam(3)+1,1,tam(3),3) ;
imgR =zeros(tam(1),tam(2)) ;
FimgR=zeros(tam(1),tam(2)) ;

% coordenadas de las ventanas de cada canal
v=coefRose(1:tam(3),2:2*tam(4)+1,tam(3),3) ;  
w(:,:,1)=real(v(:,1:2:2*tam(4)-1)); w(:,:,2)=real(v(:,2:2:2*tam(4)));%Xmin max
w(:,:,3)=imag(v(:,1:2:2*tam(4)-1)); w(:,:,4)=imag(v(:,2:2:2*tam(4)));%Ymin max
W(:,:,1)= w(:,:,2)-w(:,:,1)+1 ; W(:,:,2)= w(:,:,4)-w(:,:,3) +1 ;
centrefft=floor(tam(1:2)/2)+1 ; %centro de la fft
ampl=coefRose(tam(3)+2:2*tam(3)+1,1:tam(4),tam(3),3) ;% determinar la amplitud con la imagen gabors
if ampl==0, ampl(:,:)=1; end

%fase de los filtros:
switch param(21)
case 0, fase=1*ones(tam(3),tam(4)); %real
case 1, fase=2*ones(tam(3),tam(4)) ; %imag
case 2, fase=0*ones(tam(3),tam(4)); %complex
case 3,   fase=0*ones(tam(3),tam(4)); %complex - o, pi/2
        % en la primera escala, 0 and pi/2 son reales:
       fase(1,1)=1 ;
       if mod(tam(4),2)==0, % hay otra escala real:
         fase(1,tam(4)/2+1)=1 ;
       end
case 4, fase=0*ones(tam(3),tam(4)); fase(1,:)=1 ; %complex+ 1st scale real
case 5, fase=(max(abs(w(:,:,2)+w(:,:,1)-2*centrefft(1)),abs( w(:,:,4)+w(:,:,3)-2*centrefft(2)))<2) ; %optimized%=1 si ventana simetrica=>imagen real (0 sino)
end
fase(tam(3),:)=1; %paso bajo y alto son siempre reales



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  (2)               reconstruccion de un solo canal           %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (min(canal)>0),  %reconstruccion canales gabores o paso bajo o paso alto
   k=canal(1) ; l=canal(2) ; 
   if paramQ(6)==2, %CSF por pesos
  %         coefG(:,:,k,l)=coefG(:,:,k,l)/csf(k,l) ;%inutil si es solo un canal
   end  
   vect(1)=ceil((tam(1)+1)/2)-w(k,l,1) ; %shift para conservar el centro
   vect(2)=ceil((tam(2)+1)/2)-w(k,l,3) ;
   canalr=zeros(tam(1:2)) ; 
   canalr(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))= conj(coefRose(1:W(k,l,1),1:W(k,l,2),k,l)).*shiftFT(fft2(coefG(1:W(k,l,1),1:W(k,l,2),k,l)),vect)/ampl(k,l) ;%*4^(k*max([k,l]~=[tam(3) 2]));%HF!=LF para cumplementar
   if fase(k,l)==0,%si ==1 ventana simetrica: ya tiene los 2 gabores simetrico
        FimgR=canalr+symetriefft(canalr,1) ;
    else
        FimgR=canalr;
    end
end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       (3)          reconstruccion de todos los canales                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if canal<=0           % reconstruccion de cada canal (gabors+luminancia)+paso alto
 for k=1:tam(3)
   for l=1:tam(4) %1+(tam(4) -1)*(k<tam(3))
     if W(k,l,1)>1,   %si  hay filtro para este [k,l]
        if paramQ(6)==2, %CSF por pesos
           coefG(:,:,k,l)=coefG(:,:,k,l); %/csf(k,l) ;
        end   
        vect(1)=ceil((tam(1)+1)/2)-w(k,l,1)  ; %shift para conservar el centro
        vect(2)=ceil((tam(2)+1)/2)-w(k,l,3)  ;
        canalr=zeros(tam(1:2)) ;
        canalr(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))= conj(coefRose(1:W(k,l,1),1:W(k,l,2),k,l)).*shiftFT(fft2(coefG(1:W(k,l,1),1:W(k,l,2),k,l)),vect)/ampl(k,l) ; 
%        canal= coefRose(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4),k,l).*shiftFT(fft2(coefG(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4),k,l)),vect)*4^(k*max([k,l]~=[tam(3) 2]));%HF!=LF para cumplementar
        if fase(k,l)==0,%si ==1 ventana simetrica: ya tiene los 2 gabores simetrico
               canalr=canalr+symetriefft(canalr,1) ;
        end
        FimgR=FimgR+canalr ;
     end %W
   end %l
 end %k
end%if


%%%%%%%%%%%%%%%%%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       (4)            Modificaciones y opciones                      %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dividir por CSF y la roseta
if (size(csf)==size(FimgR))&(paramQ(6)==1),
   FimgR=FimgR./csf ; %CSF por multiplicacion del espectro
end 
if divR==1,
   roseta(roseta==0)=1 ; % si la roseta no esta definida en algunos puntos, no dividir en ellos
   FimgR =FimgR./roseta ;
end

%recortar el sobremuestreo FT img
if sizefft(1)<tam(1), 
    disp('Se recorta el sobremuestro')
  FimgR=FimgR(ceil((tam(1)+2-sizefft(1))/2:(tam(1)+sizefft(1))/2),ceil((tam(2)+2-sizefft(2))/2:(tam(2)+sizefft(2))/2)) ;
end

%reconstruccion ifft
imgR=ifft2(ifftshift(FimgR)) ;
imgRb=imgR ;

if (canal>=0)&(param(23)<3),%quitar el borde
   tb=ceil(2^(tam(3)-3))*param(23);  %tamaño del borde de cada lado de la imagen[1,1,2,4,8]
   imgR=imgR(1+tb:tam(1)-tb,1+tb:tam(2)-tb) ;
end


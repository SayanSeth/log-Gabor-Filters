function [coefRose,roseta,gabors] = logGaborFilters(Nf,No,sizefft,param)
%
% [coefRose,roseta,gabors] = logGaborFilters(Nf,No,sizefft,param)
% 
% This funtion constructs the log-Gabor transform
% designed by Sylvain Fischer (2002-2007).
%
% INPUT:
%
%   Nf - number of scales
%   No - number of orientations
%   sizefft - size of the Fourier domain. E.g. sizefft=[256 256].
%             Note that it has to be equal to the image dimensions.
%   param - contains parameters to configurate the log-Gabor transform.
%           it should be declared as a Global variable.
%
%   param(1) - Threshold for windowing filter response(e.g. =0.2).
%              It is used for subsampling sub-bands.
%              If subsampling is not of interest assign any value<0.
%
%   param(2) - First scale frequecy. E.g. param(2)=0.5.
%
%   param(3) - Useless.
%
%   param(4) - Interscale ratio. For octaves aspect set param(4)=2.
%
%   param(5:8) - Useless.
%
%   param(9) = 1 improves an overall flat response (do not tested with 
%              high-pass filters) or 0 does nothing.
%
%   param(10) - Global normalization for another improvement in
%                the overall flat response.
%             = 1 every sub-band is equally normalized,
%             = 2 high-pass filters are favored,
%             = 3 normalization is done during decomposition,
%             = 4 normalization is done during reconstruction or
%             = 0 nothing.
%
%   param(11) = 1 rotates even scales if set to 1,
%             = 2 rotates and enlarges even scales or
%             = 0 nothing.
%
%   param(12) = 1 extends domain with a border of size 2^Nf, 
%             = 0 nothing
%
%   param(13) = 3 for log-Gabor filters (other filter shape are disabled).
%
%   param(14) - Nf+1 number of scale (repeated for internal use).
%
%   param(15) - No number of orientations (repeated for internal use).
%
%   param(16) > 1  Subsampling factor if required.
%
%   param(17) - Subsampling type: 0 none, 1 Nyquist, 2 sparse, >2 sparse 2^N.
%
%   param(18) - Periodic continuity: 0 none, 1 all, 2 except 1st scale horiz/vert.
%
%   param(19) - High pass filter before first scale:
%             = 0 none,
%             = 1 isotropic (useless),
%             = 2 as sum of orientations of the highest scale or
%             = 3 filling-in. 
%   param(20) - Low pass filter:
%             = 0 none,
%             = 1 gassian,
%             = 2 isotropic (useless),
%             = 3 as sum of orientations of the lowest scale or
%             = 4 filling-in.
%   param(21) - Filter Phase:
%             = 0 real,
%             = 1 imaginary,
%             = 2 complex,
%             = 3 1st scale pi/2 real,
%             = 4 first scale real, rest is complex,
%             = 5 two first scales real or
%             = 6 reel when possible,
%
%   param(22) = 0 First scale normal,
%             = 1 built as 2 highpass filters (horizontal/vertical). 
%
%   param(23) - image border: 1 equal, 2 double, 3 none.
%
% OUTPUT:
%
%   coefRose(xf,yf,sc,or) - 4D matrix containing log-Gabor coefficients.
%             - xf, yf: fourier cartesian coordinates,
%             - or: orientation index from 1:No,
%             - sc: scale index from 1:Nf+1 (sc=Nf+1 and or=1 points the lowpass filter).
%
%   roseta - sum of all filters.
%
%   gabors (additional) - log-Gabor filters arrangement in the spatial domain (complex-valued).
%
%
% Recommended Parameters:
%       global param;
%       param = [0.2 0.5 0 2 0 0 0 0 0 1 1 0 3 Nf+1 No 0 5 2 0 3 2 0 3];
%
% Note logGaborFilters.m requires: logKernel.m, polarcoor.m and symetriefft.m
%
% For detailed aspects of the filter construction consult:
%
% S. Fischer, F. Sroubek, L. Perrinet, R. Redondo, and G. Cristóbal. “Self-invertible
% log-Gabor wavelets”. Int. Journal of Computer Vision, 75 (2), p. 231-246, 2007.
%
% This function has been upgraded from the original
% file rosette.m written by Sylvain Fischer (2002-2007).


%param(10) no esta muy claro!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      (1)     Definition de los parametros                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%penalizacionHF=1 ;%penalizar para que salgan menos coeficientes del canal HF (1 sino)

%global param
param(13) = 3;
tipo = param(13); %log-gabor por defecto

if max(sizefft)<=1, sizefft=[256 256];  end
if param(12)==1, border=2^Nf ;  %anadir un border de tamano 2^Nf
else border=0 ;
end

sizefft=sizefft+border*[1 1];% +mod(sizefft+1,2)
if param(16)>1,   %si >1 : coeficiente con el cual se va a submuestrear  la imagen antes de su proceso
  tam(1:2)= ceil([(sizefft(1)*param(16)) (sizefft(2)*param(16))]/2)*2+[1 1] ; 
else tam(1:2)=sizefft ;    % sin sobremuestreo 
end
tam(1:4)=round([tam(1), tam(2), Nf+1, No] );
if tam(4)<3,  %tam(4) tiene que estar superior a 3 para que kepa toda la informacion
    tam(4)=3 ;
end
umbralG=param(1) ;  % umbral de corte de los filtros multiresolicion (las lineas y columnas cuyas valores son todas inferiores al umbral seran puestas a cero.)
maxfrec=param(2) ;  % frecuencia de la primera escala (~0.5)
df=param(4) ;       % ratio interescala : 2 para un octavo
sigmaX=param(5) ;   % para Gabor y Vanillia : ancho de banda en X
sigmaY=param(6) ;   % lo mismo en Y
sharpness=param(7) ; % para alguno filtros : esctrecho de banda
rot= param(11) ;    % 1 si girar las escalas pares %2 si girar orientaciones pares
frec = sizefft(1)*maxfrec ;  
% param(10)         %1 si dividir los filtros por la roseta para tener una respuesta plana
divR=(param(10)==4) ;%1 si dividir los canal por la roseta a la reconstruccion (0 sino)
divD=(param(10)==3) ;%1 si dividir a la decomposicion  (0 sino)
if param(17)==0, % sin submuestreo
    param(1)=0 ; % no cortar los filtros.
end

coefRose=zeros(tam) ;
w=zeros(tam(3),tam(4),4) ;
centro=ones(tam(3),tam(4),2) ;
param(18)==1 ;

%%%%%%%%%%%%
% paramQ global variable. better to be removed before ending
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       (2)         Definicion de los filtros                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% altas frecuencias con gabores  + luminacia
escala=zeros([tam(1),tam(2),tam(3)+1]) ;  %=suma de las orientaciones de cada escala
canales=zeros(tam(1:2)) ;
%definir las coordenadas polares :
centrefft=floor(tam(1:2)/2)+1 ; %centro de la fft
x=(1:tam(1))'*ones(1,tam(2))- centrefft(1) ;
y=ones(tam(1),1)*(1:tam(2)) - centrefft(2);
rho=sqrt(x.^2 + y.^2) ; 
for k=1+(param(22)==2):tam(3)  %si se usa 2 orient en la 1era escala, la 1era escala se define al final
  do=mod(k-1,2)*mod(rot,2)*.5  ;     % rotacion de la escalas pares si rot=1 (0 si rot =0 o 2)
  frec=sizefft(1)*maxfrec/(df)^(k-1) ;
  for l=1: tam(4) 
      if (k<tam(3))&(l>No), % gabores
          break ;   % si no hay mas orientaciones, salir (por si No<3)
      end      
      if k==tam(3),         % filtro paso-bajo o alto
           roseta=sum(escala(:,:,1:tam(3)-1),3) ;
           switch l,  %tipo de filtro : PA, PB, gabor
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%     paso bajo    %%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           case 1     %metodo antiguo : a conservar para los no-logKernel   
               switch param(20) % tipo de paso-bajo %tipo
               case {0} ,
                    gab=zeros(tam(1:2)) ;
               case  1, %gaussian,Prolate centrado en 0.....%tipo{2,4,5,6,10,11,12}
                    longPB=roseta(ceil(tam(1)/2+.5),ceil(tam(2)/2+.5):tam(2)) ;
                    [tmp,longPB]=max(longPB>.5*max(max(roseta))) ;      
                    gab=canalPB(tam,frec,k,longPB,param,tipo); %*max(max(roseta)) ;
               case 2, %filtro sin orientacion {3,7,8} %solo filtros log-polar
                    gab=canal([tam(1),tam(2),tam(3),1],frec,0,k,l,param,tipo) ;     
                    centro(k,l,1:2)=ceil(tam(1:2)/2) ;
                    %rellenar el centro  :
                    longPB=gab(ceil(tam(1)/2+.5),ceil(tam(2)/2+.5):tam(2)) ;
                    [tmp,longPB]=max(longPB) ;      % da el sigmaPB
                    gab(rho<longPB)=1 ;
               case 3, %PB as an entire scale : suma de las orientaciones
                     gab=zeros(tam(1:2)) ; 
                     for k2=0:2,%suma de 3 escalas
                       do=mod(k+k2-1,2)*mod(rot,2)*.5  ; 
                       for l2=1:No,
                        if (rot==2)&(mod(l2,2)==0),  % configuracion hexagonal 2
                            frec2=frec*sqrt(df)/df^k2  ; %las orientaciones pares tienen una frecuencia aumentada.
                        else ,          frec2=frec/df^k2 ;
                        end 
                        orient=pi/No*(l2-1+do) ;
                        canal1=canal([tam(1),tam(2),tam(3),No],frec2,orient,k,l2,param,tipo) ;   
                        gab=gab+canal1.^2+symetriefft(canal1,0).^2 ;
                       end
                     end
                    gab=sqrt(gab) ;maxgab=max(max(gab)) ;
                    gab((rho<frec/sqrt(df))&(gab<maxgab*.96))=maxgab ; % relleno. A mejorar si hay rot==2
                    %se pone el umbral a .90 puesto que 2 escalas son ahnadidas, solo el centro del PB a a estar bajo el
                    %umbral, la normalizacion lo hara subrir sin lastimar los filtros de escala que tienen muy poca amplitude ahi
               case 4  %por relleno    %  case {2,4,5,6,10,11,12} %otros
                          gab=relleno(roseta,'PB')  ;
               end  %end PB
               %el paso bajo tiene que ser de tamano tam/2^(Nf-1)
               if (param(17)==1)|(param(17)==3), %downsampling 2n
                    vd2=tam(1:2)/2^(tam(3)-1);%%forzar el tamano del paso bajo:cortarlo
                    %se podria quitar el +1, pero ya no seria par? (simetrico?):
                    v=[centrefft(1)-vd2(1)+1, centrefft(1)+vd2(1)-1  centrefft(2)-vd2(2)+1 centrefft(2)+vd2(2)-1] ;%nuevas ventanas
                    gabc=zeros(tam(1:2));
                    gabc(v(1):v(2),v(3):v(4))=gab(v(1):v(2),v(3):v(4)) ; gab=gabc ;
               end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
%%%%%%%%%%%%%%%%%%    paso alto   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 2  
               switch param(19),
               case {0,4}, % sin filtro paso alto
                   gab=zeros(tam(1:2)) ;
                    break  %salir del bucle   {for l=1:tam(4)}
               case 1, % filtro sin orientacion %  switch tipo      %  case {3,7,8} %filtros logpolar
                     gab=canal([tam(1),tam(2),tam(3),1],maxfrec*sizefft(1)*df,0,k,l,param,tipo) ;       
                    %rellenar por fuera
                     long=diag(gab(ceil(tam(1)/2+.5):tam(1),ceil(tam(2)/2+.5):tam(2)))  ; 
                     [tmp,long]=max(long) ;     
                     gab(rho>long*sqrt(2))=1 ; %sqrt(2) porque es diagonal 
               case 2, %PA como suma de todas la orientaciones de una escala
                        gab=zeros(tam(1:2)) ;
                        frec=maxfrec*sizefft(1)*df ;
                        for k2=0:2,%suma de 3 escalas
                         do=mod(-k2-1,2)*mod(rot,2)*.5  ; 
                         for l2=1:No,
                             orient=pi/No*(l2-1+do) ; %orient=pi/No*(l2-1+mod(rot,2)*.5) ;
                             if (rot==2)&(mod(l2,2)==0),  % configuracion hexagonal 2
                                frec2=frec*sqrt(df)*df^k2  ; %las orientaciones pares tienen una frecuencia aumentada.
                            else ,   frec2=frec*df^k2 ;
                            end 
                            gab=gab+(canal([tam(1),tam(2),tam(3),No],frec2,orient,k,l2,param,tipo)).^2 ;   
                            gab=gab+(canal([tam(1),tam(2),tam(3),No],frec2,orient-pi,k,l2,param,tipo)).^2 ;   
                         end
                        end
                        gab=sqrt(gab) ;maxgab=max(max(gab)) ;
                        gab(rho>frec*(sqrt(df))^k2)=maxgab ; % relleno. A mejorar si rot==2.
               case 3  %por relleno    %  case {2,4,5,6,10,11,12} %otros
                        gab=relleno(roseta,'PA')  ;
               end  %end PA
            case 3
                  break  %no hay mas filtros para k=tam(3)
            end  %end switch l  : fin del PB y PA 
            if max(max(abs(gab)))>0 ,
              %  gab=gab/max(max(gab(2:tam(1),2:tam(2)))) ;%2:tam para preservar symetria
            end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%%%%%%%%%%%    filtros gabores    %%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     else       
        %modificacion frecuencia por configuracion hexagonales
        if (rot==2)&(mod(l,2)==0),  % configuracion hexagonal 2
           frec2=frec*sqrt(df)  ; %las orientaciones pares tienen una frecuencia aumentada.
        else 
           frec2=frec ;
        end 
        orient=pi/No*(l-1+do) ;%configuracion hexagonal 1
        gab=canal([tam(1),tam(2),tam(3),No],frec2,orient,k,l,param,tipo) ;        
        if (param(19)==4)&(k==1) %si juntar corners 
            switch l
            case 2 , %si juntar corners %juntar este filtro con el k==1, l==4
                gab=max(gab,canal([tam(1),tam(2),tam(3),No],frec2,orient+pi/2,k,l,param,tipo)) ;  
                gab=max(gab,canal([tam(1),tam(2),tam(3),No],frec2,orient-pi/2,k,l,param,tipo)) ;  
                %gab=gab+canal([tam(1),tam(2),tam(3),No],frec2,orient+pi/2,k,l,param,tipo) ;  
                %gab=gab+canal([tam(1),tam(2),tam(3),No],frec2,orient-pi/2,k,l,param,tipo) ;  
                 %las 3 esquinas: la 4a viene luego
            case 4,                gab=gab*0 ; %no poner nada
            end
        end
        [tmp,centro(k,l,1)]=max(max(gab,[],2)) ; 
        [tmp,centro(k,l,2)]=max(max(gab,[],1)) ;
    end
centro(tam(3),1,:)=centrefft ; %pasobajo
centro(tam(3),2,:)=[centrefft(1), 1] ;%paso alto
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   fase de los filtros de escala  : %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   if k<tam(3),
    switch param(21),
    case 0, %filtros reales
            gab2=canal([tam(1),tam(2),tam(3),No],frec2,orient-pi,k,l,param,tipo) ;  
            gab=max(gab,gab2) ; 
            fase(k,l)=1 ;
    case 1, %filtros imaginarios
        if k>1, 
            gab2=-canal([tam(1),tam(2),tam(3),No],frec2,orient-pi,k,l,param,tipo) ;  
            gab=gab+gab2 ;
        else %primera escala. hacer un shift de medio pixel en el espacio para la parte imaginaria.
            gab=canal([tam(1),tam(2),tam(3),No],frec2,orient,k,l,param,tipo,1) ;  
            gab2=canal([tam(1),tam(2),tam(3),No],frec2,orient-pi,k,l,param,tipo,1) ;  
            gab=gab-gab2 ; %filtro imaginario
        end
        fase(k,l)=2 ;
    case 2, %filtros complejos
        if k==1, %primera escala. hacer un shift de medio pixel en el espacio para la parte imaginaria.
            gab2=canal([tam(1),tam(2),tam(3),No],frec2,orient-pi,k,l,param,tipo) ;  
            gabr=max(gab,gab2)/2 ; %filtro real
            gab=canal([tam(1),tam(2),tam(3),No],frec2,orient,k,l,param,tipo,1) ;  
            gab2=canal([tam(1),tam(2),tam(3),No],frec2,orient-pi,k,l,param,tipo,1) ;  
            gabi=(gab-gab2)/2 ; %filtro imaginario
            gab=gabr+gabi ;
            %gab=gabr+gabi/1.2253 ;%proba para filtro diagonales
        end
        fase(k,l)=0 ;
    case {3,4,5}, %1eras escalas real, otras complejas.
        if (param(21)==3)&(mod(orient,pi/2)==0)&(k==1),
            fase(k,l)=1 ;   
        else %case 4,5 y 3 oblicuos   
            fase(k,l)=(k<param(21)-2) ;
        end
         if fase(k,l)==1, %reel
             gab2=canal([tam(1),tam(2),tam(3),No],frec2,orient-pi,k,l,param,tipo) ;  
             gab=max(gab,gab2) ; 
         end
    case 6, %fase optimizada para menos redundancia 
       %fase(k,l)= min(abs(v(2:2:4)+v(1:2:3)-2*centrefft)<2);  %=1 si ventana simetrica=>imagen real (0 sino)
       gab2=zeros(size(gab)) ;
       switch param(17),
       case {0,1,5},
           fase(k,l)=1 ; %real
       case {3,4},
           if k==1, fase(k,l)=1;
           else fase(k,l)=0;
           end
       case 2,  fase(k,l)=0;%complex
       end
       if fase(k,l)==1 , 
           %ventana simetrica: coger los 2 gabores simetricos->senal real (no para LP ni HP)
            gab2=canal([tam(1),tam(2),tam(3),No],frec2,orient-pi,k,l,param,tipo) ;  
            %  figure(11) ;surf(gab) ;figure(12) ; surf(gab2) ;          
       end  %sino ventana asimetrica: solo 1 gabor 
       gab=max(gab,gab2) ;   %    gab=gab+gab2 ;
    end%switch
  else %k==tam(3)
      fase(k,l)=1 ; %los paso-alto y bajo son reales;
  end %if
      
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   tipo de submuestreo  : %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     [a,b]=find(abs(gab)>umbralG) ;          % cortar la columnas y filas que esten debajo de un umbral
     v=[min(a) max(a) min(b) max(b)] ;   % determinar la ventana del filtro
     if max(max(abs(gab)))>0 %si hay un filtro :
        switch param(17) %tipo de downsampling
        case 0 %no downsampling
            v=[1 tam(1) 1 tam(2)] ; %sin submuestreo
        case 1 % Nyquist downsampling
            vd2=2.^(ceil(log2([max(abs(v(1:2)-centrefft(1)))   max(abs(v(3:4)-centrefft(2))) ]))) ;% tamano de ventanas
            v=[centrefft(1)-vd2(1), centrefft(1)+vd2(1)-1  centrefft(2)-vd2(2) centrefft(2)+vd2(2)-1] ;%nuevas ventanas
%      case 2 %sparse downsampling
        case {3,4,5} % ^2 sparse downsampling: para imagenes de 2^Mx2^N pixeles
            vd=[v(2)-v(1)+1, v(4)-v(3)+1 ] ;  %tamano actual de ventana
            switch param(17),
            case 3, vd2= 2.^ceil([log2(v(2)-v(1))  log2(v(4)-v(3))]) ; % tamano de ventana 2^n posible
            case 4, 
                vd2= sizefft(1:2)/2^(k-1); %.^2
                if [k,l]==[tam(3) 1], vd2=vd2*2 ;%paso bajo              
                end
            case 5, vd2= tam(1:2)/max(1,2^(k-2));
                if [k,l]==[tam(3) 1], vd2=vd2*2 ;%paso bajo              
                end
            end
            if [k,l]==[tam(3) 2], vd2=sizefft(1:2) ;%paso alto:ventana entera
            end
            vdd=((vd2-vd)/2) ; % demi diferencia de tamano de ventana.
            v=[v(1)-ceil(vdd(1)), v(2)+floor(vdd(1)), v(3)-ceil(vdd(2)), v(4)+floor(vdd(2))] ; %nueva ventana
            if min(tam(1:2)-vd2)==0, %una de las dimensiones de ventana es tan grande que la imagen
                v=[1 tam(1) 1 tam(2)]; %coger la fft entera y se aprovechara de la symetria.
            end
            if v(1)<1,               v(2)=v(2)-v(1)+1 ; v(1)=1 ;            end %desplazar si se queda fuera
            if v(2)>sizefft(1),     v(1)=v(1)-v(2)+sizefft(1) ; v(2)=sizefft(1) ;            end
            if v(3)<1,                 v(4)=v(4)-v(3)+1 ; v(3)=1 ;            end
            if v(4)>sizefft(2),     v(3)=v(3)-v(4)+sizefft(2) ; v(4)=sizefft(2) ;            end
            %ajustar mas finamente la ventana:
          if k<tam(3), %menos para el paso bajo
            if v(2)<tam(1),
              %  while max(abs(gab(v(1),:)),[],2)<max(abs(gab(v(2)+1,:)),[],2),
                while sum(abs(gab(v(1),:)),2)<sum(abs(gab(v(2)+1,:)),2),
                    v(1:2)=v(1:2)+1 ;
                    if v(2)==tam(1), break; end
                end
            end
            if v(1)>1,
                while sum(abs(gab(v(1)-1,:)),2)>sum(abs(gab(v(2),:)),2),
                    v(1:2)=v(1:2)-1 ;
                    if v(1)==1, break; end
                end
            end    
            if v(4)<tam(2),
                while sum(abs(gab(:,v(3))),1)<sum(abs(gab(:,v(4)+1)),1),
                    v(3:4)=v(3:4)+1 ;
                    if v(4)==tam(2), break; end
                end
            end
            if v(3)>1,
                while sum(abs(gab(:,v(3)-1)),1)>sum(abs(gab(:,v(4))),1),
                    v(3:4)=v(3:4)-1 ;
                    if v(3)==1, break; end
                end
            end    
          end%menos para el paso bajo
        end
       % reel(k,l)= min(abs(v(2:2:4)+v(1:2:3)-2*centrefft)<2);  %=1 si ventana simetrica=>imagen real (0 sino)
             
        gabc=zeros(tam(1:2));     
        gabc(v(1):v(2),v(3):v(4))=gab(v(1):v(2),v(3):v(4)); %gabc = filtro cortado para downsampling
        if [k,l]==[tam(3),1], %paso bajo
            gabc(v(1),:)=0 ; gabc(:,v(3))=0 ; %para la symetria quitar las primeras lineas y columnas del paso bajo
        end
           

        
%%%%%% %%%%%%%     guardar el filtro (k,l) : %%%%%%%%%%%%%%%%%%
        if (param(22)<2)|(k<tam(3))|(l~=2), % 1er escala normal.
        w(k,l,:)=v ;  
        coefRose(1:v(2)-v(1)+1,1:v(4)-v(3)+1,k,l)=gabc(v(1):v(2),v(3):v(4)) ; %guardar el filtro  
        coefRose(k,2*l:2*l+1,Nf+1,3)=[v(1)+i*v(3)  v(2)+i*v(4)] ;  %guardara las coordenadas de ventana del filtro (min max) %%%a hacer fuera del bucle???        
        else % 1er escala con 2 orientaciones.% crea 2 filtros de 1era escala a partir del paso-alto
           w(1,:,:)=0 ; w(1,1,:)=v ; w(1,2,:)=v ; 
           [rho,teta]=polarcoor(gabc) ;
           gabcH=sqrt(.5*cos(2*teta)+.5).*gabc ;
           gabcV=sqrt(-.5*cos(2*teta)+.5).*gabc ;
           coefRose(1:v(2)-v(1)+1,1:v(4)-v(3)+1,1,1)=gabcH(v(1):v(2),v(3):v(4)) ; %guardar el filtro  
           coefRose(1,2:3,Nf+1,3)=[v(1)+i*v(3)  v(2)+i*v(4)] ;  %guardara las coordenadas de ventana del filtro (min max)     
           coefRose(1:v(2)-v(1)+1,1:v(4)-v(3)+1,1,2)=gabcV(v(1):v(2),v(3):v(4)) ; %guardar el filtro  
           coefRose(1,4:5,Nf+1,3)=[v(1)+i*v(3)  v(2)+i*v(4)] ;  %guardara las coordenadas de ventana del filtro (min max)           
        end% 1er escala con 2 orientaciones.
 %%%%%%%%%%%%%% calculo de la roseta y escalas: a hacer en una funccion a parte. %%%%%%%%%%%%%%%%%%%%%
        if fase(k,l)>0,            invgab=zeros(tam(1:2)) ;    
        else  invgab=symetriefft(gabc,1);% filtro symetrico para complementar 
        end
        if [k l]==[tam(3) 2]  % paso alto
                escala(:,:,k+1) = gabc.*conj(gabc) + invgab.*conj(invgab) ;
        else% gabores orientables y paso bajo
                escala(:,:,k) = escala(:,:,k) + gabc.*conj(gabc) + invgab.*conj(invgab) ; 
        end
     end %if hay un filtro
  end  % end for l
end %end for k
roseta=sum(escala,3) ; %%%%%%%%%%%%%%%??????????????????????
if (param(19)==0)|(param(19)==4) % no hay paso alto por separado
    rosetaLF=sum(escala(:,:,2:tam(3)),3) ; %todo menos las altas freq.
    rosetaHF=escala(:,:,1) ; %todo menos las altas freq.
else
    rosetaLF=sum(escala(:,:,1:tam(3)),3) ; %todo menos las altas freq.
    rosetaHF=escala(:,:,tam(3)+1) ; %todo menos las altas freq.
end
 W(:,:,1)= w(:,:,2)-w(:,:,1)+1; W(:,:,2)= w(:,:,4)-w(:,:,3) +1 ; %ventanas de los filtros
% reel=(max(abs(w(:,:,2)+w(:,:,1)-2*centrefft(1)),abs(w(:,:,4)+w(:,:,3)-2*centrefft(2)))<2) %ya calculado

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      (3)      Modification de los filtros para mejorar la roseta       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Opcion Param(9) : poner coeficientes a las escalas para aplanar la roseta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param(9)==1,  %no comprabado con un paso alto (param(19)>0)
    tm=tam(3)+(param(19)>0)-(param(20)==0) ;
    cEsc=ones(tm,1) ; maxEsc=zeros(tm,1) ; m=0; 
    while max(abs(maxEsc-1))>0.02 & m<8,
        m=m+1 ;
        roseta=sum(escala,3) ;
        for k=1:tm,    % solo para roseta redonda   %+1 si hay un paso-alto
            if (k==tm)&(param(19)>0),
                maxEsc(k)=roseta(centro(tam(3),2,1),centro(tam(3),2,2)) ; 
            else
                maxEsc(k)=roseta(centro(k,1,1),centro(k,1,2)) ; 
            end
            maxEsc(1)=1; %no bajar la altura de las HF.
            escala(:,:,k)=escala(:,:,k)/(maxEsc(k)^.8) ;
            cEsc(k)=cEsc(k)/(maxEsc(k)^.4 ) ;   %escala ya esta al cuadrado
        end
        disp(maxEsc )  
    end
    for k=1:tam(3)-1
        coefRose(:,:,k,:)=coefRose(:,:,k,:)*cEsc(k) ;
    end
    coefRose(:,:,tam(3),1)=coefRose(:,:,tam(3),1)*cEsc(tam(3)) ;
    roseta=sum(escala,3) ;   %roseta valida para la reconstrucion. YA NO
    if (param(19)==0)|(param(19)==4) % no hay paso alto por separado
        rosetaLF=sum(escala(:,:,2:tam(3)),3) ; %todo menos las altas freq.
        rosetaHF=escala(:,:,1) ; % altas freq.   
    else
        rosetaLF=sum(escala(:,:,1:tam(3)),3) ; %todo menos las altas freq.
        rosetaHF=escala(:,:,tam(3)+1) ; %altas freq.
    end
end

% opcion param(10) : modificar los canales para una roseta plana 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param(10)==1, %rebajar cada filtro igual
  roseta1=roseta;
  roseta1(abs(roseta1)<1e-19)=1e-19 ;
  roseta1=sqrt(roseta1) ;
  roseta=zeros(tam(1:2)) ;
%  figure(15) ; showGray(roseta1) ;
  for k=1:tam(3)  %dividir los filtros por la (roseta1)^.5 y calcular la nueva roseta%no se cambia las HF
    for l=1: tam(4), %1 + (No-1)*(k<tam(3)) + (param(19)>0)*(k==tam(3)) 
       if W(k,l,1)>1,   %si  hay filtro para este [k,l]
  %     coefRose(1:W(k,l,1),1:W(k,l,2),k,l)=coefRose(1:W(k,l,1),1:W(k,l,2),  %     k,l)./((roseta1(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))).^.5) ;
       coefRose(1:W(k,l,1),1:W(k,l,2),k,l)=coefRose(1:W(k,l,1),1:W(k,l,2),k,l)./roseta1(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4)) ;
       canal2=zeros(tam(1:2)) ;
       canal2(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))=coefRose(1:W(k,l,1),1:W(k,l,2),k,l) ;
       roseta=roseta+canal2.*conj(canal2) ;
       if fase(k,l)==0, %filtro complejo, complementar por el symetrico
           symcanal=symetriefft(canal2,1);
           roseta=roseta+symcanal.*conj(symcanal) ;
       end %if
%       canales(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))=coefRose(1:W(k,l,1),1:W(k,l,2),k,l) ; % imagen de la superosition de los filtros frecuenciales
       end
    end
  end    
end

if param(10)==2, %poner ventaja al paso alto o 1st scale: no rebajarle
  roseta1=roseta;
  roseta1(abs(roseta1)<1e-19)=1e-19 ;
  rosetaLF(abs(rosetaLF)<1e-19)=1e-19 ;
  rosetaHF(abs(rosetaHF)<1e-19)=1e-19 ;
  %aumentar el nivel del Paso alto hasta rellenar roseta a 1.
  rosetaUp=max(min(1,roseta1),rosetaHF) ; %inferior a 1:->crecimiento de los filtros HF  %>1->disminucion
  if (param(19)==0)|(param(19)==4) % si no hay paso alto por separado
     HPtype=2 ;
     k=1 ;     rosetaHF =zeros(tam(1:2)) ;
     for l=1: tam(4), %high pass por la primera escala.
          if W(k,l,1)>1,   %si  hay filtro para este [k,l]
          coefRose(1:W(k,l,1),1:W(k,l,2),k,l)=coefRose(1:W(k,l,1),1:W(k,l,2),k,l)./((rosetaUp(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))).^.5) ; 
          canal2=zeros(tam(1:2)) ;
          canal2(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))=coefRose(1:W(k,l,1),1:W(k,l,2),k,l).^2 ;
          rosetaHF=rosetaHF+canal2 ;
%          rosetaHF=rosetaHF+coefRose(:,:,k,l).^2 ; %todo las altas freq.   
          if max([tam(1)-W(k,l,1),tam(2)-W(k,l,2)])>0, %si hay submuestreo: complementar
              canal2=symetriefft(canal2,1) ;
              rosetaHF=rosetaHF+canal2 ;     
          end
          end
     end%end for l
  else
     HPtype=1 ;%high pass en el canal tam(3),2
     k=tam(3) ; l=2 ; %el paso-alto%se supone que no esta submuestreado
     coefRose(1:W(k,l,1),1:W(k,l,2),k,l)=coefRose(1:W(k,l,1),1:W(k,l,2),k,l)./((rosetaUp(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))).^.5) ;  
     rosetaHF=coefRose(:,:,tam(3),2).*conj(coefRose(:,:,tam(3),2)) ; %todo menos las altas freq.   
  end%if no hay canal HF separado
  roseta=zeros(tam(1:2)) ;
  for k=HPtype:tam(3)  %dividir los filtros por la (roseta1)^.5 y calcular la nueva roseta%no se cambia las HF
     for l=1: tam(4), %1 + (No-1)*(k<tam(3)) + (param(19)>0)*(k==tam(3)) 
         if (k==tam(3))&(l>1),      break; end %no se disminuye las HF
         if W(k,l,1)>1,   %si  hay filtro para este [k,l]
  %     coefRose(1:W(k,l,1),1:W(k,l,2),k,l)=coefRose(1:W(k,l,1),1:W(k,l,2),  %     k,l)./((roseta1(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))).^.5) ;
       correction=sqrt(abs(1-rosetaHF)./rosetaLF) ;
       coefRose(1:W(k,l,1),1:W(k,l,2),k,l)=coefRose(1:W(k,l,1),1:W(k,l,2),k,l).*correction(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4)) ;
       canal2=zeros(tam(1:2)) ;
       canal2(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))=coefRose(1:W(k,l,1),1:W(k,l,2),k,l) ;
       roseta=roseta+canal2.*conj(canal2) ;
       %roseta(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))=roseta(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))+coefRose(1:W(k,l,1),1:W(k,l,2),k,l).^2 ;
       if fase(k,l)==0, %complementar por el symetrico
           symcanal=symetriefft(canal2,1);
           roseta=roseta+symcanal.*conj(symcanal) ;
           %roseta=roseta+symetriefft(canal2,1).^2 ;
       end %if
    %   canales(w(k,l,1):w(k,l,2),w(k,l,3):w(k,l,4))=coefRose(1:W(k,l,1),1:W(k,l,2),k,l) ; % imagen de la superosition de los filtros frecuenciales
        end%if W>1
    end%for l
  end%for k
 roseta=roseta+rosetaHF ;
end%

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    (4)   Registrar parametros y crear imagenes  de los filtros espaciales       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     guardar los tamanos de la fft, border y tam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%coefRose(tam(3)+1,1,tam(3),3)=border ; % guardar el tamano del border
coefRose(Nf+2,2,Nf+1,3)=tam(1)+i*tam(2) ; %guardar tam
coefRose(Nf+2,4:5,Nf+1,3)=sizefft ; %guardar tamano de la imagen+1+border
%tammax(1)=  max([max(max(W(:,:,1))),2*tam(4)+2,17])  ;   % reducir el tamano de coefRose quitando las lineas y
%tammax(2)=  max([max(max(W(:,:,2))),2*tam(4)+2,17])  ;   %   columnas nulas
%coefRose= coefRose(1:tammax(1),1:tammax(2),:,:) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     gabors = imagen de los filtros en el dominio espacial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im=zeros(sizefft) ;
csf=ones(sizefft);
[coefL,csf] = dcmp(im,csf,coefRose,divD,roseta) ; %sin CSF
% calcular la amplitud de cada gabor en la imagen reconstruida y guardarlo en coefrose:
ampl= zeros(tam(3:4)) ;
for k=1:Nf+1,  %gabores
  for l=1:No,
      if W(k,l,1)>1,
        coefG=coefL ;  coefG(ceil(W(k,l,1)/2),ceil(W(k,l,2)/2),k,l)=1;%1
        gab=rcmp(coefG,divR,csf,roseta,coefRose,0) ;
        ampl(k,l)=max(max(max(real(gab)))-min(min(real(gab))),max(max(imag(gab)))-min(min(imag(gab))));
      end
  end
end
if (1==0)&(W(Nf+1,2,1)<2), %si no hay paso alto, calcular la amplitude de la 1era escala de manera diferente
    amplsc1=ampl(1,:) ;
    amplsc1(amplsc1==0)=inf ; %para evitar tomar en cuenta filtros inexistantes
    ampl(1,:)=min(min(ampl(1,:))) ; 
    ampl(1,:)=min(amplsc1).*(ampl(1,:)>0) ; %para la primera escala se promedia las amplitudes porque no hay sufficiemmente precision para el calculo
end
%ampl(2,:)=min(ampl(2,:)) ; %a quitar
%guardar las amplitudes calculadas
ampl;
coefRose(tam(3)+2:2*tam(3)+1,1:tam(4),tam(3),3)=ampl(:,:) ; 
%imagen final de gabores teniendo en cuenta las amplitudes
%poner coeficientes no nulos en cada canal :
[coefG,csf] = dcmp(im,csf,coefRose,divD,roseta) ; %sin CSF
for k=1:Nf,
    for l=1:No
        coefG(ceil(W(k,l,1)/(Nf+1)*(Nf-k+1.5)),ceil(W(k,l,2)/No*(l-.5)),k,l)=1;%1
    end
end
coefG(ceil(W(tam(3),1,1)/(Nf+1)*.5),ceil(W(tam(3),1,2)*.25),Nf+1,1)=1;%.2; %en el paso-bajo : a 1/4
coefG(ceil(W(tam(3),2,1)/(Nf+1)*.5),ceil(W(tam(3),2,2)*.75),Nf+1,2)=1;%1.5; %paso-alto
%coefG(:,:,tam(3),3)=coefRose(:,:,tam(3),3) ; showPyr2(coefG) ; figure ;
param23sav=param(23) ; param(23)=3 ; %para las imagenes no quitar el borde
gabors=real(rcmp(coefG,divR,ones(sizefft),roseta,coefRose,0)) ;
gabors=gabors+i*real(rcmp(i*coefG,divR,ones(sizefft),roseta,coefRose,0)) ;
param(23)=param23sav ;
clear paramQ
%min(min(gabors)) , max(max(gabors))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      calcular la forma de las maskaras :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[coefMask] = dcmp(gabors,ones(sizefft),coefRose,divD,roseta) ;
%showPyr(coefMask,-4) ;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       (5)         funcciones de definicion de los fitros          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                   Paso banda                     %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gab=canal(tam,frec,orient,k,l,param,tipo,spaceshift)
% gab=canal(tam,frec,orient,k,l,param,tipo)
if exist('spaceshift')==0, spaceshift=0;
end
  
%diag=(mod(orient,pi/2)>.1)&(mod(orient,pi/2)<pi/2-.1)  ;%no va bien si se usa >19 orientaciones
diag=(mod(orient,pi/2)>.03)&(mod(orient,pi/2)<pi/2-.03)  ;
contiP=(k<3) & (((param(18)==1) & ((k>1)|(diag)|(param(21)~=3))) | ((param(18)==2)&((diag)|(k>1))));
faseshift=tam(1) ;
if contiP %continuidad periodica
    tam(1:2)=tam(1:2)*3 ; 
end
switch tipo,
case 2, gab=gabor(tam(1)+1,orient,frec,param(5),param(6)*tam(4)/4) ;
case {3,7,8}, gab=logKernel(tam(1),orient,frec,	tam(4),tipo,0) ;  
case 4, gab=chocoPie(tam(1)+1,orient,frec*4/tam(1),.5*tam(4),0.5) ;
case 5 ,gab=vanilla(tam(1)+1,orient,frec,param(5),param(6)*tam(4)/4) ;
case 6, %gab=stretch(tam(1)+1,orient,frec,param(7),0) ; %now steerable.
    imgs=zeros(tam(1:2)) ; imgs(1,1)=1 ;
    [PYR, INDICES, STEERMTX, HARMONICS] = buildFullSFpyr2(imgs,tam(3),tam(4)-1, 1) ; %faire la transf. HP
    l2=l+tam(4)*mod(k-1,2) ;
    k2=ceil(k/2) ; %disp([k2 l2]) ;
    [gab] = spyrBand(PYR,INDICES,k2,l2) ; %band pass 
    gab2=real(fftshift(fft2(-i*gab))) ; 
    tmg=size(gab2) ;%tamano del steerable
    gab=zeros(tam(1:2)) ;
    gab(floor(tam(1)/2-tmg(1)/2)+1:floor(tam(1)/2+tmg(1)/2),floor(tam(2)/2-tmg(2)/2)+1:floor(tam(2)/2+tmg(2)/2))=gab2 ; % enlever le sous echantillonnage
case {10,11,12}, gab=paving(tam,l,	k,tipo) ;
end
 % showGray(gab, 11) ; r=input('?','s') ;
gab=gab(1:tam(1),1:tam(2)); %disminuir de una linea y columna en su caso

%disp([k,l,frec,orient/pi]) ; figure(11) ; showGray(gab) ;

if (spaceshift==1), %&(k==1), %desplazar de medio pixel en el espacio.
    x=([1:tam(1)]'-floor(tam(1)/2)-1)*ones(1,tam(2)) ;
    y=ones(tam(1),1)*([1:tam(2)]-floor(tam(2)/2)-1) ;
    orx=-sin(mod(orient,pi)) ;%modulo pi: independiente de si es +/- pi.
    ory=cos(mod(orient,pi)) ;
    faseshift=faseshift*1/max(abs(orx),abs(ory)) ; %maxima frequencia en la direccion del Gabor.
    faseshift=exp(i*pi*(x*orx+y*ory)/faseshift) ;%solo para direccion real    
    gab=gab.*faseshift ;   
    
end
%figure(12) ; showGray(imag(faseshift)) ;

%gab=gab/max(max(gab(2:tam(1),2:tam(2)))) ;%2:tam para preservar symetria                %           end 
if contiP,
    tam(1:2)=tam(1:2)/3 ; 
    gab1=zeros(tam(1:2)) ; 
%    gabplie=zeros(tam(1:2)) ; %seulement pour les graphiques
    for kx=1:3, %repliegar
        switch kx,
        case 1, attx= (tam(1)+1-[1:tam(1)]')*ones(1,tam(2))/tam(1) ;%Â¡Â¡attencion a la asymetria de la fft!!
        case 2, attx=zeros(tam(1:2)) ;
        case 3, attx= ([1:tam(1)]'-1)*ones(1,tam(2))/tam(1) ;
        end
        for ky=1:3
            switch ky,
            case 1, atty= ones(tam(1),1)*(tam(2)+1-[1:tam(2)])/tam(2) ;
            case 2, atty=zeros(tam(1:2)) ;
            case 3, atty= ones(tam(1),1)*([1:tam(2)]-1)/tam(2) ;
            end
            att=pi*3*sqrt(attx.^2+atty.^2); %vale pi a un  tercio de la imagen
            att(att>pi)=pi;
                        %con atenuacion proporcional al alejarse
          %  gabplie=gabplie+gab(1+(kx-1)*tam(1):kx*tam(1),1+(ky-1)*tam(2):ky*tam(2)) ;%seulement pour les graphiques
            gab1=gab1+gab(1+(kx-1)*tam(1):kx*tam(1),1+(ky-1)*tam(2):ky*tam(2)).*((cos(att)+1)/2).^(1) ;
            %gab1=max(gab1,gab(1+(kx-1)*tam(1):kx*tam(1),1+(ky-1)*tam(2):ky*tam(2))) ;
        end
    end
    
%graphique:
if (1==0)&(k==1)&(l==2)&(spaceshift==0),
    figure ; showGray(gab) ; contour(gab,.6) ;
   figure ; showGray(gabplie) ;contour(gabplie,.6) ;
    figure ; showGray(gab1) ;contour(gab1,.6) ; 
    figure ;
end
    
   % gab1=gab(1+tam(1):2*tam(1),1+tam(2):2*tam(2)) ;
    gab=gab1 ; 
  %  figure(10) ; showGray(gab) ; r=input('next?','s') ;
  
  
   % figure(13) ; showGray(imag(gab)) ;
    %gabb=gab/max(max(abs(gab))) ; 
    %dessin(:,:,1)=real(gab) ;dessin(:,:,2)=real(gab) ;
    %dessin(:,:,3)=abs(imag(gab)) ;
    %figure(14) ; showGray(dessin) ;
%    figure(14) ; showGray(real(faseshift)) ;
%    figure(15) ; showGray(imag(faseshift)) ;
% r=input('next?','s') ;   
end
% figure(12) ; showGray(real(gab)) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                   Paso bajo                     %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gab=canalPB(tam,frec,k,longPB,param,tipo)
 centrefft=floor(tam(1:2)/2)+1 ; %centro de la fft
X=((1:tam(1)) - centrefft(1)).^2 ;
Y=((1:tam(2)) - centrefft(2)).^2 ;
switch tipo
case {2},      % gaussiana    
  sigmaPB = longPB/sqrt(log(2)) ; %tam(1)*param(2)/pow2(tam(3)-2)*param(3)*2 ; 
  gab= exp(-(ones(tam(1),1)*Y + X'*ones(1,tam(2)))/sigmaPB^2 ) ;
case {4,5,6}  % prolate
    milong =floor(tam(1)/2) ;           % probar que anchura tiene el filtro
    val=zeros(tam(1),1) ;
    E=dpss(2*milong+1,2) ; val(1:milong+1)=E(milong+1:2*milong+1,1) ;
    [tmp,longmedia]=max(val<sqrt(.5)*max(val)) ;       
    milong=round(milong*longPB/longmedia) ;  % ajustar con la anchura requerida
    val=zeros(tam(1),1) ;               % y recalcular el filtro
    E=dpss(2*milong+1,2) ; val(1:milong+1)=E(milong+1:2*milong+1,1) ;   
    gab=val(round(sqrt(ones(tam(1),1)*Y + X'*ones(1,tam(2))))+1) ; 
case {3,7,8}  % log Kernel
    gab=logKernel(tam(1),0,frec,tam(4),tipo,1) ; 
case {10,11,12}, gab=paving(tam,0,tam(3),tipo) ; 
end  
gab =gab/max(max(gab)) ;  % sin normalizacion,  max gabor= 1



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                   relleno                       %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filtro=relleno(roseta,AB) 
%rellenar las partes no cubiertas de la roseta en 
%un paso alto si AB='PA'
%o un paso bajo  si AB='PB' 

[rho,teta,teta2,X,Y]=polarcoor(roseta) ;
tam=size(roseta) ;
% centrefft=floor(tam(1:2)/2)+1 ; %centro de la fft
x=(1:tam(1))'*ones(1,tam(2)) ;
y=ones(tam(1),1)*(1:tam(2)) ;
%X=x - centrefft(1);  %tam(1)/2  - 1/2 ;
%Y=y - centrefft(2); %tam(2)/2  - 1/2 ;
%rho=(Y.^2 + X.^2).^.5 ; 
rho(rho==0)=.00001 ;

filtro=1-roseta ;
filtro(filtro<0)=0 ;
filtro=sqrt(filtro) ;

dX=round(X./rho);
dY=round(Y./rho);
x3=x-dX;
y3=y-dY ;
dX=dX.*(x<tam(1)).*(x>1) ;
dY=dY.*(y<tam(2)).*(y>1) ;
x2=x+dX ; 
y2=y+dY ; 
nelim=1;
switch AB
case 'PA' %coger solo la parte creciente, poner a cero la parte decreciente
    filtro(rho<4)=0 ; % quitar el centro
    while nelim>0
            eliminar=filtro(c1d(x,y,tam))>filtro(c1d(x2,y2,tam)) ;  
            nelim=sum(sum(eliminar)) ;
        
            filtro(eliminar==1)=0 ;
    end
case 'PB' %poner a cero la parte creciente
      while nelim>0
            eliminar=filtro(c1d(x,y,tam))>filtro(c1d(x3,y3,tam)) ;  
            nelim=sum(sum(eliminar)) ;
            filtro(eliminar==1)=0 ;   
            disp(nelim) 
    end
end

function coor1d=c1d(x,y,tam) ;
coor1d=x+tam(1)*(y-1) ;

%function img=linetoimg(L,tam) ;
%for k=1:tam(2),
%    img(1:tam(1),k)=L((k-1)*tam(1)+1:k*tam(1)) ;
%end





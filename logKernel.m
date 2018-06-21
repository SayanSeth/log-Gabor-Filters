function logK = logKernel(tam,teta0,f0,No,tipo,PB) 
% logK = logKernel(tam,teta0,f0,No,tipo,PB) ;
%construye un filtro en coordenadas log-polares
%tam: dimensiones de la imagen :have to bee even
% teta0: orientacion preferiada del filtro
% f0: frequencia preferida
% No: numero de orientaciones
% tipo: tipo de filtros a construir
% PB= 1 si es el paso bajo o 0 sino
% Ajustar el argumento del prolate directamente en el programa

global param ;

%teta0=mod(teta0+pi,2*pi)-pi ;%partenece a [-pi, pi[
teta0=mod(teta0+pi/2,2*pi)-pi/2 ;%partenece a [-pi/2,3/2*pi[
%teta0=mod(teta0,2*pi) ;%partenece a [0,2*pi[

logK=zeros(tam,tam) ;
c=floor(tam/2)+1;  %centro de la fft %c= ceil((tam + 1)/2) ; 
if No>1,
    Dteta=pi/No ;
else  % si hay una sola orientacion, el ancho angular del filtro es infinito
    Dteta=1e14 ;
end

[f,teta,teta2]=polarcoor([tam tam]) ;
%X=(1:tam)-c ;
%f=sqrt((X.^2)'*ones(1,tam) + ones(tam,1)*(X.^2)) ;
f(c,c) = 0.00001 ;  %triche
%teta=-asin((X'*ones(1,tam))./f) ;
%teta2=teta+2*pi;
%teta2(:,1:c-1)=pi-teta(:,1:c-1) ;
%teta(:,1:c-1)=-pi-teta(:,1:c-1) ;
   
   
%   teta0
%   figure(2) ;mesh(teta') ;xlabel('x') ;ylabel('y') ;
%   figure(3) ;mesh(teta2') ;xlabel('x') ;ylabel('y') ;
%   figure(4) ; mesh(mod(abs(teta-teta0),2*pi)') ;xlabel('x') ;ylabel('y') ;
%  figure(5) ; mesh(mod(abs(teta2-teta0),2*pi)') ; xlabel('x') ; ylabel('y') ;
%  r=input('next?') ;
tet=min(mod(abs(teta-teta0),2*pi),mod(abs(teta2-teta0),2*pi))/Dteta ;
%tet=mod(abs(teta-teta0),2*pi)/Dteta ;
%tet=(teta-teta0)/Dteta ;
rho=log2(f/f0) ;

% calculo de la distancia al centro del filtro a construir para cada punto del plano de Fourier
if PB==0 , 
    switch param(11),   % si shift o hexagonal
    case 0 ,  dist=(rho.^2 +tet.^2).^.5;
    case {1} ,  rho=3^.5/2*rho ; % roseta shift ;
        if param(8)==0 ,  dist=(rho.^2 +tet.^2).^.5;
        else % hexagonal distance en roseta shift
            beta=param(8) ;   %1/log2(2/3^.5) ~=4.82
            dist= (abs(3^.5/2*rho+.5*tet).^beta+abs(-3^.5/2*rho+.5*tet).^beta+abs(tet).^beta).^(1/beta);
        end
    case 2
            tet=3^.5/2*tet ; % paar scale shift
            if param(8)==0 ,  dist=(rho.^2 +tet.^2).^.5;
            else % hexagonal distance en roseta shift
                beta=param(8) ;   %1/log2(2/3^.5) ~=4.82
                dist= (abs(3^.5/2*tet+.5*rho).^beta+abs(-3^.5/2*tet+.5*rho).^beta+abs(rho).^beta).^(1/beta);
           end
    end
else ,    %filtro paso bajo
   switch param(11),   % si rotation entre escalas: cambia la distancia
   case {0} , dist=rho  ;  
   case {1,3} , f00=f0*2*2^(-.75^(-.5)) ;
            dist=sqrt(.75)*log2(f/f00)  ;
   case 2   %   param(11)==2 (escalas paares desplacadas)
           b=4 ; %forma x^b+y^b del baso-bajo : b=2 para forma circular
           correccionf0=2^(1/4) ; % correccion debida a que las orientationes no tienen las mismas escalas
            fPB=((abs(X)'*ones(1,tam)).^b+(ones(tam,1)*abs(X)).^b).^(1/b) ; %habia puesto sqrt
            fPB(c,c) = 0.00001 ;  %triche
            dist=log2(fPB/f0/correccionf0) ;
   end
   dist(dist<0)=0 ;
end

% dibujar el filtro basandose en la distancia definida
switch tipo
case 3,        % log Gabor:
        if param(11)==0, K=2*log(2) ;
        else K=.996; %K=1.056; %experimentalemente 06agosto2004: minimiza la diferencia absoluta entre max y min de la roseta (.6100334, 3.5729% )% 
            %K=.996: experimentalmente 06agosto2004: minimiza el porcentage de diferencia entre max y min de la rosetan (abs dif=6.27511, 3.47752%)%
            %K=1; %  (max-min)/min (roseta en LF)= 1.27%  28fev05
            %K=1.2431070547 ; %en roseta shift si %  3X^4-X^3-X+.5=0  y  K=-2logX pero no da mejor resultado
                end 
   %     K=2.3 ;
    logK=exp(-K*dist.^2);
    % if PB==0, % correcion de amplituda si no paso bajo por el solapamiento de los filtros vecinos
    %   logK=logK/(1  + 6* exp(-2*K)) ;
    % end
case 6,     %steerable
    logK=(.5+.5*cos(rho/2*pi)).*(.5+.5*cos(tet/2*pi)) ;
    logK(abs(rho)>2)=0 ;
    logK(abs(tet)>2)=0 ; % no es asi, es cos^3 (voir simoncelli 1991) 
    
case 7,        % log Prolate
 % dist2(dist2>1.5)=1.5 ;
 % VER SI SE PUEDE PONER MAS ESTRECHO
  c=513 ;   X=zeros(2*c-1) ;
  E=dpss(c,1.5) ;  %2 : lobulo secundario-50dB, 1.5:-40dB, 1:-30 dB
  X(ceil(c/2):floor(c*1.5))=E(:,1)/max(max(E(:,1))) ;
  %encontrar el X=2^(-.5)  :la suma de 2 filtros vale 1 en su medio (en dist=.5)
  %m,X2]=max(X(c+1:2*c-1)<2^(-.5))  ;
  %logK=X(round(dist*2*X2+c)) ;  
  % la suma de 3 filtros vale 1 entre los 3 : en dist = sqrt(3)/3=.5774
  %poner el factore entre 1.0 (dominio de transform. pequeno ) y 1.1 (solapamiento uniforme de los filtros)
    [m,X3]=max(X(c+1:2*c-1)<1/sqrt(3)*1.06) ;
    logK=X(round(dist*sqrt(3)*X3+c)) ;  
    
case 8,        % external 
        global filtroext ; 
        tamfiltext=max(size(filtroext)) ; 
        tamfiltU= round(tamfiltext/4) ;  % CUIDADO con cambiar el valor
        %  [tmp,xmin]=max(filtroext>0) ;
        %  [tmp,xmax]=max((filtroext+((1:tamfiltext)<xmin))==0) ; 
        xmid= (1+tamfiltext)/2 ;
         xb=tamfiltext/2+tamfiltU/2+.5 ; 
        dist(dist>1.4)=1.4 ;
        logK=(filtroext(round(xmid+dist*2*(xb-xmid)))) ; 
        % suavizar :
        suma=conv2(abs(logK)>0,ones(3,3),'same') ;
        suma(suma==0)=-1;
        for n=1:3,
            nonezero=(logK~=0) ;
            logK=conv2(logK,ones(3,3),'same');
            logK=logK.*nonezero./suma ;
        end
end
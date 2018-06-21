function [rho,teta,teta2,X,Y]=polarcoor(tam) ;
% [rho,teta,teta2,X,Y]=polarcoor(tam) ;
% define las coordenadas polares
% teta in [-3pi/2,pi/2]
% teta2 in [pi/2, 3pi/2]
 
if length(tam)>4, tam=size(tam) ;
end
centrefft=floor(tam(1:2)/2)+1 ; %centro de la fft
x=(1:tam(1))'*ones(1,tam(2)) ;
y=ones(tam(1),1)*(1:tam(2)) ;
X=x - centrefft(1);  %tam(1)/2  - 1/2 ;
Y=y - centrefft(2); %tam(2)/2  - 1/2 ;
rho=(Y.^2 + X.^2).^.5 ; 
f=rho ;
f(rho==0)=.00001 ;

teta=-asin(X./f) ;
teta2=teta+2*pi;
teta2(:,1:centrefft(2)-1)=pi-teta(:,1:centrefft(2)-1) ;
teta(:,1:centrefft(2)-1)=-pi-teta(:,1:centrefft(2)-1) ;


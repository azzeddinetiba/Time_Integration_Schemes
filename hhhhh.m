clear ;
clc
E=200e9;
S=0.002;
L=2;
m=100;
T=0.08;
t=linspace(0,T,800);
x=2;
pho=5000 ;
p=-3500;

ustat=(p*x/(E*S));
SOL=zeros(1,size(t,2));


for i=1:size(t,2)
    h=0;
for n=1:m 

w=@(n)(2*n-1) * pi/(2*L) * sqrt( E /pho );
q=@(n)(2*n-1) * pi/(2*L) ;
ui=@(n,x,T)(((-1)^(n-1))/(2*n-1)^2)*cos(w(n)*T)*sin(q(n)*x) ;
h=h+ui(n,x,t(i)) ;
end
SOL(1,i)=ustat*0+(8*p*L/((pi^2)*E*S))*h;
end

%SOL=ustat-(8*p*L/((pi^2)*E*S))*SOL ;

plot(t,SOL);
%fprintf(' %d ',U);


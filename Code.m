clear ;
clc

T=10^(-6);
%T=0.000001;
FORCE=3500;
pho=5000;
S=0.002;
E=200e9;
L=2;
M=pho*S*L;


figure('name','freq propres')
theor=zeros(1,13);
for i=1:13
    theor(i)=(i*(pi/L)*sqrt(E/pho))/(2*pi);
end

plot(linspace(1,13,13),theor)
hold on

for incr=0:5
    
Nmesh=15+7*incr;

MAT=zeros(Nmesh,Nmesh);
K=zeros(Nmesh,Nmesh);
F=zeros(Nmesh,1);

me=M/(Nmesh-1);
mat_elem_M=(me/6)*[2 1;1 2];
le=L/(Nmesh-1);
mat_elem_K=((E*S)/le)*[1 -1;-1 1];

for k=1:Nmesh-1
for i=1:2
for j=1:2
I=k+i-1;
J=k+j-1;
K(I,J)=K(I,J)+ mat_elem_K(i,j);
end
end
end

for k=1:Nmesh-1
for i=1:2
for j=1:2
I=k+i-1;
J=k+j-1;
MAT(I,J)=MAT(I,J)+ mat_elem_M(i,j);
end
end
end

for k=1:Nmesh-1
elemFi = SM_elem;
for i=1:2
I=k+i-1;
F(I)=F(I)+ elemFi(i);
end
end

F(1:Nmesh-1)=zeros(Nmesh-1,1);

for k=1:Nmesh
MAT(k,1)=0;
MAT(1,k)=0;
end
MAT(1,1)=1;
for k=1:Nmesh
K(k,1)=0;
K(1,k)=0;
end
K(1,1)=1;


Nt=floor((0.08)/T);

A=(T^2)*K+MAT;
A2=eye(Nmesh)+A\MAT;
U2=(T^2)*(A2\(A\F));

SOL=zeros(Nmesh,Nt);
SOL(:,2)=U2;



 for j=3:Nt
U=A\(((T^2)*F)-(MAT*SOL(:,j-2))+(2*MAT*SOL(:,j-1)));
SOL(:,j)=U;
 end

 
%AFF=zeros(Nt,1);
%for j=1:Nt
%AFF(j)=SOL(Nmesh,j);
%end
 
%figure('name', 'Graphe de déplacement de x=L');
%X=linspace(0,0.08,Nt);
%plot(X,AFF)

K=K(2:end-1,2:end-1); MAT=MAT(2:end-1,2:end-1);
[X,Lam] = eig(K,MAT);

Lam = Lam.^(1/2);
Lam=Lam./(2*pi);



mef=zeros(1,13);
for i=1:13
    mef(i)=Lam(i,i);
end

plot(linspace(1,13,13),mef)
hold on


%figure('name', 'Deplacement de la poutre');
    %colormap hsv;
    
    %for i=2:Nt
    %for i=2:5
    %DEPL=SOL(:,i);
    %Z=[DEPL';DEPL';DEPL';DEPL']';Y=linspace(0,0.1,4);
    %surf(Y,linspace(0,L,15),Z,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    %colorbar;
    %daspect([2 1*abs(DEPL(15)*10^5) 100*(L)]); axis off; title( sprintf('t = %f s', i*T) );
    %pause(1.4)
    %end
    
end

legend('theor','15','22','29','36','43','50')
hold off;




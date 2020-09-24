%Differences finies [0 ; 0.008s]

clear ;
clc

MAT=zeros(15,15);
K=zeros(15,15);
F=zeros(15,1);
T=10^(-5);
%T=0.000001;
FORCE=3500;
pho=5000;
S=0.002;
E=200e9;
L=2;
M=pho*S*L;
me=M/14;
mat_elem_M=(me/6)*[2 1;1 2];
le=L/14;
mat_elem_K=((E*S)/le)*[1 -1;-1 1];

for k=1:14
for i=1:2
for j=1:2
I=k+i-1;
J=k+j-1;
K(I,J)=K(I,J)+ mat_elem_K(i,j);
end
end
end

for k=1:14
for i=1:2
for j=1:2
I=k+i-1;
J=k+j-1;
MAT(I,J)=MAT(I,J)+ mat_elem_M(i,j);
end
end
end

for k=1:14
elemFi = SM_elem;
for i=1:2
I=k+i-1;
F(I)=F(I)+ elemFi(i);
end
end

F(1:14)=zeros(14,1);

for k=1:15
MAT(k,1)=0;
MAT(1,k)=0;
end
MAT(1,1)=1;
for k=1:15
K(k,1)=0;
K(1,k)=0;
end
K(1,1)=1;


Nt=floor((0.008)/T);

A=(T^2)*K+MAT;
A2=eye(15)+A\MAT;
U2=(T^2)*(A2\(A\F));

SOL=zeros(15,Nt);
SOL(:,2)=U2;



 for j=3:Nt
U=A\(((T^2)*F)-(MAT*SOL(:,j-2))+(2*MAT*SOL(:,j-1)));
SOL(:,j)=U;
 end

 
AFF=zeros(Nt,1);
for j=1:Nt
AFF(j)=SOL(15,j);
end
 
figure('name', 'Graphe de déplacement de x=L');
X=linspace(0,0.008,Nt);
plot(X,AFF)

figure('name', 'Deplacement de la poutre');
    colormap hsv;
    
    %for i=2:Nt
    for i=2:5
    DEPL=SOL(:,i);
    Z=[DEPL';DEPL';DEPL';DEPL']';Y=linspace(0,0.1,4);
    surf(Y,linspace(0,L,15),Z,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    colorbar;
    daspect([2 1*abs(DEPL(15)*10^5) 100*(L)]); axis off; title( sprintf('t = %f s', i*T) );
    pause(1.4)
    end
    

hold off;



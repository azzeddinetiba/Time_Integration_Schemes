%  Schéma explicite Diff Finies centrées   %

clear ;
clc

MAT=zeros(15,15);
K=zeros(15,15);
F=zeros(15,1);
%T=0.00001310251;
T=10^(-5)+0.000003;
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


Nt=floor((0.08)/T);

%L'initiation de la matrice Solution
SOL=zeros(15,Nt);

h=T^2;

for i=3:Nt
X=(MAT/h)\(((2/h)*MAT-K)*SOL(:,i-1)-((1/h)*MAT)*SOL(:,i-2)+F);
SOL(:,i)=X;
end


AFF=zeros(Nt,1);
for j=1:Nt
AFF(j)=SOL(15,j);
end
 
figure('name', 'Graphe de déplacement de x=L');
X=linspace(0,0.08,Nt);
plot(X,AFF) 

%AFF=zeros(Nt,1);
%for j=1:Nt
%AFF(j)=SOL(30,j);
%end
%figure('name', 'Graphe de vitesse à x=L');
%plot(X,AFF) 

%AFF=zeros(Nt,1);
%for j=1:Nt
%AFF(j)=SOL(3*15,j);
%end
%figure('name', 'Graphe de l accélération à x=L');
%plot(X,AFF)


    



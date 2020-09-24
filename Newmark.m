function SOL=Newmark(NELT,XAFF)

MAT=zeros(NELT+1,NELT+1);
K=zeros(NELT+1,NELT+1);
F=zeros(NELT+1,1);
T=0.0001;
FORCE=3500;
pho=5000;
S=0.002;
E=200e9;
L=2;
M=pho*S*L;
me=M/NELT;
mat_elem_M=(me/6)*[2 1;1 2];
le=L/NELT;
mat_elem_K=((E*S)/le)*[1 -1;-1 1];

for k=1:NELT
for i=1:2
for j=1:2
I=k+i-1;
J=k+j-1;
K(I,J)=K(I,J)+ mat_elem_K(i,j);
end
end
end

for k=1:NELT
for i=1:2
for j=1:2
I=k+i-1;
J=k+j-1;
MAT(I,J)=MAT(I,J)+ mat_elem_M(i,j);
end
end
end

for k=1:NELT
elemFi = SM_elem;
for i=1:2
I=k+i-1;
F(I)=F(I)+ elemFi(i);
end
end

F(1:NELT)=zeros(NELT,1);

for k=1:NELT+1
MAT(k,1)=0;
MAT(1,k)=0;
end
MAT(1,1)=1;
for k=1:NELT+1
K(k,1)=0;
K(1,k)=0;
end
K(1,1)=1;


Nt=(0.08)/0.0001;

SOL=zeros(3*(NELT+1),Nt);
Beta=0.7;


for i=2:Nt
A=SOL(1:NELT+1,i-1)+T*SOL(NELT+1+1:(NELT+1)*2,i-1)+(1-Beta)*0.5*(T^2)*SOL((NELT+1)*2+1:3*(NELT+1),i-1);
B=SOL(NELT+1+1:(NELT+1)*2,i-1)+(1-Beta)*T*SOL((NELT+1)*2+1:3*(NELT+1),i-1);
Xderivsecond=-(MAT+0.5*Beta*(T^2)*K)\(K*A-F);
SOL((NELT+1)*2+1:3*(NELT+1),i)=Xderivsecond;
X=A+0.5*Beta*(T^2)*Xderivsecond;
SOL(1:(NELT+1),i)=X;
Xderiv=B+Beta*T*Xderivsecond;
SOL(NELT+1+1:(NELT+1)*2,i)=Xderiv;
end



AFF=zeros(Nt,1);
for j=1:Nt
AFF(j)=SOL(XAFF,j);
end
 
X=linspace(0,0.08,Nt);
plot(X,AFF) 

 


end
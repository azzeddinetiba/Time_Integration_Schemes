function SOL=dynamics(NELT,XAFF)

Noeuds=NELT+1;

MAT=zeros(Noeuds,Noeuds);
K=zeros(Noeuds,Noeuds);
F=zeros(Noeuds,1);
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
    elemFi = SM_elem;
for i=1:2
for j=1:2
I=k+i-1;
J=k+j-1;
K(I,J)=K(I,J)+ mat_elem_K(i,j);
MAT(I,J)=MAT(I,J)+ mat_elem_M(i,j);
end
F(I)=F(I)+ elemFi(i);
end
end


F(1:NELT)=zeros(NELT,1);

for k=1:Noeuds
MAT(k,1)=0;
MAT(1,k)=0;
end
MAT(1,1)=1;
for k=1:Noeuds
K(k,1)=0;
K(1,k)=0;
end
K(1,1)=1;


Nt=(0.08)/0.0001;

A=(T^2)*K+MAT;
A2=eye(Noeuds)+A\MAT;
U2=(T^2)*(A2\(A\F));

SOL=zeros(Noeuds,Nt);
SOL(:,2)=U2;

 for j=3:Nt
U=A\(((T^2)*F)-(MAT*SOL(:,j-2))+(2*MAT*SOL(:,j-1)));
SOL(:,j)=U;
 end



 
AFF=zeros(Nt,1);
for j=1:Nt
AFF(j)=SOL(XAFF,j);
end
 
X=linspace(0,0.08,Nt);
plot(X,AFF)

end
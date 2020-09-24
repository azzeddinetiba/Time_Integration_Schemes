function  y = phi1_P1(x, x1, x2)
%calcule la fonction de forme phi1 telle que: 
%     phi1(x1)=1 et phi1(x2)=0
%               Ti
%       |---------------|
%       x1              x2
% phi1(x)=a*x+b et ab=[a;b]
ab=[x1 1;x2 1]\[1;0];
    y = polyval(ab,x);
end
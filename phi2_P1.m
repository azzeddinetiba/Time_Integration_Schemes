function  y = phi2_P1(x, x1, x2)
%calcule la fonction de forme phi2 telle que: 
%     phi2(x1)=0 et phi2(x2)=1
%               Ti
%       |---------------|
%       x1              x2
% phi1(x)=a*x+b et ab=[a;b]
ab=[x1 1;x2 1]\[0;1];
    y = polyval(ab,x);
end
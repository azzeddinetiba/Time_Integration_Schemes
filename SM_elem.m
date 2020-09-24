function [SM] = SM_elem

%xm=(k*le+(k-1)*le)/2;
%SM=[-3500*(le/6)*(phi1_P1((k-1)*le,(k-1)*le,k*le)+phi1_P1(k*le,(k-1)*le,k*le)+4*phi1_P1(xm,(k-1)*le,k*le));
 %   -3500*(le/6)*(phi2_P1((k-1)*le,(k-1)*le,k*le)+phi2_P1(k*le,(k-1)*le,k*le)+4*phi2_P1(xm,(k-1)*le,k*le))];

 SM=[-3500;-3500];
 
end


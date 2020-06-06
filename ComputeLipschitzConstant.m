function [LipschitzConst] = ComputeLipschitzConstant(Mesh,PhiComponent,...
    FL2Norm,SourceNumUpdate)
%   Compute L = \int_{Q} A^2(Q) dQ, where A(Q) = \phi(Q)/G_{P_i}; 

Integral = zeros(SourceNumUpdate,1);
D = sin(Mesh.ThetaQ);
dx = Mesh.ThetaQLine';
dy = Mesh.PsiQLine';
for i = 1:SourceNumUpdate
    Integrand = ((PhiComponent(:,:,i)./FL2Norm(i)).^4).*D;
    Integral(i) = trapz(dy,trapz(dx,Integrand,2));
end
LipschitzConst = sum(trapz(dy,trapz(dx,Integrand,2)));

end


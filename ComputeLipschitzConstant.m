function [LipschitzConst] = ComputeLipschitzConstant(Mesh,PhiComponent,...
    FL2Norm,SourceNumUpdate)
%   Compute L = \int_{Q} A^2(Q) dQ, where A(Q) = \phi(Q)/G_{P_i}; 

D = sin(Mesh.ThetaQ);
dx = Mesh.ThetaQLine;
dy = Mesh.PsiQLine;
Ai = zeros(size(PhiComponent));
for i = 1:SourceNumUpdate
    Ai(:,:,i) = PhiComponent(:,:,i);
end
AiSum = sum(Ai,3);
AiNorm = sqrt(sum(Ai.^2,3));
Integrand = AiNorm.*AiSum.*D;
LipschitzConst = abs(trapz(dy,trapz(dx,Integrand,2)));

end


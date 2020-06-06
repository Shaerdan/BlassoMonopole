function [FL2Norm] = L2NormF(Locations,Mesh)

%%   Compute the L2norm of the forward integral kernel F(Q,S):

%%   Initialization:
SourceNum = length(Locations)/3;
SourceRadius = Locations(1:SourceNum);
SourceTheta = Locations(SourceNum+1:2*SourceNum);
SourcePsi = Locations(2*SourceNum+1:3*SourceNum);
dx = Mesh.ThetaQLine';
dy = Mesh.PsiQLine';
FL2Norm   = zeros(1,SourceNum);

%%   Compute the kernel F and its L2norm:
for i=1:SourceNum
    CosGamma = cos(SourceTheta(i)).*cos(Mesh.ThetaQ)+sin(SourceTheta(i)).*...
    sin(Mesh.ThetaQ).*cos(SourcePsi(i)-Mesh.PsiQ);
    DistancePQ = sqrt(1 + SourceRadius(i).^2 - 2* SourceRadius(i).*CosGamma);
    F = (2./DistancePQ) - log(1-SourceRadius(i).*CosGamma + DistancePQ);
    Integrand = (F.^2).*sin(Mesh.ThetaQ);
    FL2Norm(i) = sqrt(trapz(dy, trapz(dx,Integrand,2)));
end
end


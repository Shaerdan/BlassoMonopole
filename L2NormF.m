function [FL2Norm] = L2NormF(Locations,Mesh)

%%   Compute the L2norm of the forward integral kernel F(Q,S):

%%   Initialization:
global LogSwitch
SourceNum = length(Locations)/3;
SourceRadius = Locations(1:SourceNum);
SourceTheta = Locations(SourceNum+1:2*SourceNum);
SourcePsi = Locations(2*SourceNum+1:3*SourceNum);
FL2Norm   = zeros(1,SourceNum);

%%   Compute the kernel F and its L2norm:
for i=1:SourceNum
    CosGamma = cos(SourceTheta(i))*cos(Mesh.ThetaQ(:))+sin(SourceTheta(i))*...
    sin(Mesh.ThetaQ(:)).*cos(SourcePsi(i)-Mesh.PsiQ(:));
    DistancePQ = sqrt(1 + SourceRadius(i)^2 - 2*SourceRadius(i)*CosGamma);
    F = (2./DistancePQ) - LogSwitch*log(1-SourceRadius(i)*CosGamma + DistancePQ);
    FL2Norm(i) = norm(F,2);
end
end


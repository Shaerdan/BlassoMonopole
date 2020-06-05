function [Eta] = ObjEta(X,Lambda,Mesh,P)
%%   Compute the objective function -(1/2)*Eta;
%   Eta -- ref Equation (45), Page 6 of the report: 
%          Eta = Eta(Mesh.S) = integral (
%          F_Normalized(Mesh.Q,Mesh.S)*P(Mesh.Q)dMesh.Q ),
%          dMesh.Q = sin(Mesh.ThetaQ) * dMesh.ThetaQ * dMesh.PsiQ. 
%   F  --- The integral kernel, ref Equation (45), Page 6 of the report.
%   FL2Norm --- L2Norm of F;
%% 

dx = Mesh.ThetaLine';
dy = Mesh.PsiLine';
CosGamma = cos(X(2)).*cos(Mesh.ThetaQ)+sin(X(2)).*...
    sin(Mesh.ThetaQ).*cos(X(3)- Mesh.PsiQ);
DistanceSQ  = sqrt(1 + X(1).^2 - 2* X(1).*CosGamma);
F      = (2./DistanceSQ) - log(1-X(1).*CosGamma + DistanceSQ);
Integrand0 = (F.^2)*sin(Mesh.ThetaQ);
FL2Norm = sqrt(trapz(Mesh.PsiLine', trapz(Mesh.ThetaLine',Integrand0,2)));
F_Normalized  = F./FL2Norm;
Integrand1 = F_Normalized.*P.*sin(Mesh.ThetaQ);
Eta = abs(trapz(dy, trapz(dx,Integrand1,2)))/Lambda;
Eta = -0.5*Eta^2;

end


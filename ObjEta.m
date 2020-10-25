function [Eta] = ObjEta(X,Lambda,Mesh,P,SwitchNorm)
%%   Compute the objective function -(1/2)*Eta;
%   Eta -- ref Equation (45), Page 6 of the report: 
%          Eta = Eta(Mesh.S) = integral (
%          F_Normalized(Mesh.Q,Mesh.S)*P(Mesh.Q)dMesh.Q ),
%          dMesh.Q = sin(Mesh.ThetaQ) * dMesh.ThetaQ * dMesh.PsiQ. 
%   F  --- The integral kernel, ref Equation (45), Page 6 of the report.
%   FL2Norm --- L2Norm of F;
%% 

dx = Mesh.ThetaQLine;
dy = Mesh.PsiQLine;
CosGamma = cos(X(2)).*cos(Mesh.ThetaQ)+sin(X(2)).*...
    sin(Mesh.ThetaQ).*cos(X(3)- Mesh.PsiQ);
DistanceSQ  = sqrt(1 + X(1)^2 - 2* X(1)*CosGamma);
F      = (2./DistanceSQ) - log(1-X(1)*CosGamma + DistanceSQ);
Integrand0 = (F.^2).*sin(Mesh.ThetaQ);
if (SwitchNorm == 1)
    FL2Norm = sqrt(trapz(dy, trapz(dx,Integrand0,2)));
else
    FL2Norm = 1;
end
F_Normalized  = F./FL2Norm;
Integrand1 = F_Normalized.*P.*sin(Mesh.ThetaQ);
Eta = abs(trapz(dy, trapz(dx,Integrand1,2)))/Lambda;
Eta = -0.5*Eta^2;

end


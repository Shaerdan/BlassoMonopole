function [PhiP] = Phi_P(Mesh,P)
%   Given P(x), compute Phi(P)
%%
%   PhiP -- Phi(P).
%   F    -- Integral Kernel.
%%
dx1 = Mesh.ThetaS'; dy1 = Mesh.PsiS'; dz1 = Mesh.RadiusS';
n1 = length(Mesh.ThetaQLine);
n2 = length(Mesh.PsiQLine);
SIN1 = sin(Mesh.STheta);
R1 = Mesh.SRadius.^2;
for i=1:n1
    for j=1:n2
        CosGamma = cos(Mesh.STheta).*cos(Mesh.ThetaQLine(i))+sin(Mesh.STheta).*...
            sin(Mesh.ThetaQLine(i)).*cos(Mesh.SPsi- Mesh.PsiQLine(j));
        DistanceSQ  = sqrt(1 + Mesh.SRadius.^2 - 2* Mesh.SRadius.*CosGamma);
        F      = (2./DistanceSQ) - log(1-Mesh.SRadius.*CosGamma + DistanceSQ);
        %             if (SwitchNorm == 1)
        %                 Integrand0 = (F.^2).*sin(Mesh.ThetaQ);
        %                 FL2Norm = sqrt(trapz(dy, trapz(dx,Integrand0,2)));
        %             else
        %                 FL2Norm = 1;
        %             end
        %             F_Normalized  = F./FL2Norm;
        %             Integrand1 = F_Normalized.*Q.*sin(Mesh.ThetaQ);
        Integrand0 = F.*P.*R1.*SIN1;
        PhiP(i,j) = trapz(dx1,trapz(dy1,trapz(dz1,Integrand0,3),2),1);
    end
end
end


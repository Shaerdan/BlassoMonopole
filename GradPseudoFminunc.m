function [GradPseudo] = GradPseudoFminunc(X,Measurement,Mesh,Lambda,SwitchNorm)

%%  Compute pseudo gradient of ||Estimates(X) - Measurement||_2^2 + ...
%  Lambda* ||X||_1;
%  GradPseudo = d||Estimates - Measurement||_2^2/dX + ...
%  SubGradient of Lambda*||X||_1;
%  I --- Source intensity;

%% Initialize:
SourceNum = length(X)/4;
I = X(1:SourceNum);
Locations = X(SourceNum+1:end);
RadiusPi = X(SourceNum +1:2*SourceNum);
ThetaPi = X(2*SourceNum +1:3*SourceNum);
PsiPi = X(3*SourceNum +1:4*SourceNum);
SIN1 = sin(Mesh.ThetaQ);
COS1 = cos(Mesh.ThetaQ);
dx = Mesh.ThetaQLine;
dy = Mesh.PsiQLine;
GradPseudo = zeros(1,SourceNum);
GradSub = GradPseudo;
GradDescrepency  = GradPseudo;
[PhiPi,CosGamma,DistPQ] = ComputePotentialComponent(Locations,Mesh);
if (SwitchNorm == 1)
    FL2Norm = L2NormF(Locations,Mesh);
    Estimation=PhiPi*(I./FL2Norm')';
    Descrepency = Estimation - Measurement;
    %% Compute the Pseudo gradient of the objective function:
    for i = 1:SourceNum
        
        GPi = FL2Norm(i);
        %% dObjectiveFunction/dIntensity
        GradDescrepency  = (PhiPi./GPi)'*Descrepency;
        GradSub = Lambda*sign(X);
        GradPseudo(i) =   GradDescrepency(i) + GradSub(i);
        
        %% dObjectiveFunction/dRadius
        PDotQ = RadiusPi(i)*CosGamma(:,:,i);
        Di_DistPQ = (RadiusPi(i)-CosGamma(:,:,i))./DistPQ(:,:,i);
        Di_PDotQ  = CosGamma(:,:,i);
        Di_PhiPi  = (-2./(DistPQ(:,:,i).^2)).*Di_DistPQ ...
            - (-Di_PDotQ + Di_DistPQ)./(1- PDotQ + DistPQ(:,:,i));
        Di_GPi    = trapz(dy,trapz(dx,PhiPi(:,:,i).*Di_PhiPi.*SIN1,2))/GPi;
        Di_Estimate = I(i)*(GPi*Di_PhiPi - PhiPi(:,:,i)*Di_GPi)/(GPi)^2;
        Integrand = Descrepency.*Di_Estimate.*SIN1;
        GradDescrepency(i+SourceNum)  = trapz(dy,trapz(dx,Integrand,2));
        GradPseudo(i+SourceNum) = GradDescrepency(i+SourceNum) ;
        
        %% dObjectiveFunction/dTheta
        Di_CosGamma = -sin(ThetaPi(i))*COS1+cos(ThetaPi(i))*SIN1...
            .*cos(Mesh.PsiQ - PsiPi(i));
        Di_DistPQ = -(RadiusPi(i)*Di_CosGamma./DistPQ(:,:,i));
        Di_PDotQ  = RadiusPi(i)*Di_CosGamma;
        Di_PhiPi  = (-2./(DistPQ(:,:,i).^2)).*Di_DistPQ ...
            - (-Di_PDotQ + Di_DistPQ)./(1- PDotQ + DistPQ(:,:,i));
        Di_GPi    = (trapz(dy,trapz(dx,PhiPi(:,:,i).*Di_PhiPi.*SIN1,2))/GPi);
        Di_Estimate = I(i)*(GPi*Di_PhiPi - PhiPi(:,:,i)*Di_GPi)/(GPi)^2;
        Integrand = Descrepency.*Di_Estimate.*SIN1;
        GradDescrepency(i+2*SourceNum)  = trapz(dy,trapz(dx,Integrand,2));
        GradPseudo(i+2*SourceNum) = GradDescrepency(i+2*SourceNum) ;
        
        %% dObjectiveFunction/dPsi
        Di_CosGamma =  -sin(ThetaPi(i)).*SIN1...
            .*sin(PsiPi(i) - Mesh.PsiQ);
        Di_DistPQ = -(RadiusPi(i)*Di_CosGamma./DistPQ(:,:,i));
        Di_PDotQ  = RadiusPi(i)*Di_CosGamma;
        Di_PhiPi  = (-2./(DistPQ(:,:,i).^2)).*Di_DistPQ ...
            - (-Di_PDotQ + Di_DistPQ)./(1- PDotQ + DistPQ(:,:,i));
        Di_GPi    = (trapz(dy,trapz(dx,PhiPi(:,:,i).*Di_PhiPi.*SIN1,2))/GPi);
        Di_Estimate = I(i)*(GPi*Di_PhiPi - PhiPi(:,:,i)*Di_GPi)/(GPi)^2;
        Integrand = Descrepency.*Di_Estimate.*SIN1;
        GradDescrepency(i+3*SourceNum)  = trapz(dy,trapz(dx,Integrand,2));
        GradPseudo(i+3*SourceNum) = GradDescrepency(i+3*SourceNum) ;
    end
else
    FL2Norm = 1;
    Estimation=sum(bsxfun(@times,PhiPi,reshape(I./FL2Norm',1,1,SourceNum)),3);
    Descrepency = Estimation - Measurement;
    for i = 1:SourceNum
        
        %% dObjectiveFunction/dIntensity
        Integrand_Di_f_des = (Descrepency.*(PhiPi(:,:,i)).*SIN1);
        GradDescrepency(i)  = trapz(dy,trapz(dx,Integrand_Di_f_des,2));
        GradSub(i) = Lambda*sign(X(i));
        GradPseudo(i) =   GradDescrepency(i) + GradSub(i);
        
        %% dObjectiveFunction/dRadius
        PDotQ = RadiusPi(i)*CosGamma(:,:,i);
        Di_DistPQ = (RadiusPi(i)-CosGamma(:,:,i))./DistPQ(:,:,i);
        Di_PDotQ  = CosGamma(:,:,i);
        Di_PhiPi  = (-2./(DistPQ(:,:,i).^2)).*Di_DistPQ ...
            - (-Di_PDotQ + Di_DistPQ)./(1- PDotQ + DistPQ(:,:,i));
        Di_Estimate = I(i)*(Di_PhiPi);
        Integrand = Descrepency.*Di_Estimate.*SIN1;
        GradDescrepency(i+SourceNum)  = trapz(dy,trapz(dx,Integrand,2));
        GradPseudo(i+SourceNum) = GradDescrepency(i+SourceNum) ;
        
        %% dObjectiveFunction/dTheta
        Di_CosGamma = -sin(ThetaPi(i))*COS1+cos(ThetaPi(i))*SIN1...
            .*cos(Mesh.PsiQ - PsiPi(i));
        Di_DistPQ = -(RadiusPi(i)*Di_CosGamma./DistPQ(:,:,i));
        Di_PDotQ  = RadiusPi(i)*Di_CosGamma;
        Di_PhiPi  = (-2./(DistPQ(:,:,i).^2)).*Di_DistPQ ...
            - (-Di_PDotQ + Di_DistPQ)./(1- PDotQ + DistPQ(:,:,i));
        Di_Estimate = I(i)*(Di_PhiPi);
        Integrand = Descrepency.*Di_Estimate.*SIN1;
        GradDescrepency(i+2*SourceNum)  = trapz(dy,trapz(dx,Integrand,2));
        GradPseudo(i+2*SourceNum) = GradDescrepency(i+2*SourceNum) ;
        
        %% dObjectiveFunction/dPsi
        Di_CosGamma =  -sin(ThetaPi(i)).*SIN1...
            .*sin(PsiPi(i) - Mesh.PsiQ);
        Di_DistPQ = -(RadiusPi(i)*Di_CosGamma./DistPQ(:,:,i));
        Di_PDotQ  = RadiusPi(i)*Di_CosGamma;
        Di_PhiPi  = (-2./(DistPQ(:,:,i).^2)).*Di_DistPQ ...
            - (-Di_PDotQ + Di_DistPQ)./(1- PDotQ + DistPQ(:,:,i));
        Di_Estimate = I(i)*(Di_PhiPi);
        Integrand = Descrepency.*Di_Estimate.*SIN1;
        GradDescrepency(i+3*SourceNum)  = trapz(dy,trapz(dx,Integrand,2));
        GradPseudo(i+3*SourceNum) = GradDescrepency(i+3*SourceNum) ;
    end
end

GradPseudo = GradPseudo';
end

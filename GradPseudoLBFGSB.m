function [GradPseudo] = GradPseudoLBFGSB(X,Descrepency,PhiPi,CosGamma,DistPQ,Mesh,FL2Norm, Lambda)

%%  Compute pseudo gradient of ||Estimates(X) - Measurement||_2^2 + ...
%  Lambda* ||X||_1;
%  GradPseudo = d||Estimates - Measurement||_2^2/dX + ...
%  SubGradient of Lambda*||X||_1;
%  I --- Source intensity;

%% Initialize:
SourceNum = length(X)/4;
I = X(1:SourceNum);
Locations = X(SourceNum +1:end);
A1 = sin(Mesh.ThetaQ);
A2 = cos(Mesh.ThetaQ);
dx = Mesh.ThetaQ(1,:)';
dy = Mesh.PsiQ(:,1);
GradPseudo = zeros(1,SourceNum);
GradSub = GradPseudo;
GradDescrepency  = GradPseudo;

%% Compute the Pseudo gradient of the objective function:
for i = 1:SourceNum
    %% dObjectiveFunction/dIntensity
    Integrand_Di_f_des = (Descrepency.*(PhiPi(:,:,i)/FL2Norm(i)).*A1);
    GradDescrepency(i)  = trapz(dy,trapz(dx,Integrand_Di_f_des,2));
    GradSub(i) = Lambda*sign(X(i));
    GradPseudo(i) =   GradDescrepency(i) + GradSub(i);
    
    %% dObjectiveFunction/dRadius
    PDotQ = Locations(i)*CosGamma(:,:,i);
    Di_DistPQ = (Locations(i)-CosGamma(:,:,i))./DistPQ(:,:,i);
    Di_PDotQ  = CosGamma(:,:,i);
    Di_PhiPi  = (-2./(DistPQ(:,:,i).^2)).*Di_DistPQ ...
        - (-Di_PDotQ + Di_DistPQ)./(1-PDotQ + DistPQ(:,:,i));
    Di_GPi    = (trapz(dy,trapz(dx,PhiPi(:,:,i).*Di_PhiPi,2))/FL2Norm(i));
    Di_Estimate = I(i)*(FL2Norm(i)*Di_PhiPi + PhiPi(:,:,i).*Di_GPi)/(FL2Norm(i))^2;
    Integrand = Descrepency.*Di_Estimate.*A1;
    GradDescrepency(i+SourceNum)  = trapz(dy,trapz(dx,Integrand,2));
    GradPseudo(i+SourceNum) = GradDescrepency(i+SourceNum) ;
    
    %% dObjectiveFunction/dTheta
    Di_CosGamma = -sin(Locations(i+SourceNum))*A2+cos(Locations(i+SourceNum))*A1...
        .*cos(Mesh.PsiQ - Locations(i+2*SourceNum));
    Di_DistPQ = -(Locations(i)*Di_CosGamma./DistPQ(:,:,i));
    Di_PDotQ  = Locations(i)*Di_CosGamma;
    Di_PhiPi  = (-2./(DistPQ(:,:,i).^2)).*Di_DistPQ ...
        - (-Di_PDotQ + Di_DistPQ)./(1-PDotQ + DistPQ(:,:,i));
    Di_GPi    = (trapz(dy,trapz(dx,PhiPi(:,:,i).*Di_PhiPi,2))/FL2Norm(i));
    Di_Estimate = I(i)*(FL2Norm(i)*Di_PhiPi + PhiPi(:,:,i).*Di_GPi)/(FL2Norm(i))^2;
    Integrand = Descrepency.*Di_Estimate.*A1;
    GradDescrepency(i+2*SourceNum)  = trapz(dy,trapz(dx,Integrand,2));
    GradPseudo(i+2*SourceNum) = GradDescrepency(i+2*SourceNum) ;
    
    %% dObjectiveFunction/dPsi
    Di_CosGamma =  sin(Locations(i+SourceNum)).*A1...
        .*sin(Mesh.PsiQ - Locations(i+2*SourceNum));
    Di_DistPQ = -(Locations(i)*Di_CosGamma./DistPQ(:,:,i));
    Di_PDotQ  = Locations(i)*Di_CosGamma;
    Di_PhiPi  = (-2./(DistPQ(:,:,i).^2)).*Di_DistPQ ...
        - (-Di_PDotQ + Di_DistPQ)./(1-PDotQ + DistPQ(:,:,i));
    Di_GPi    = (trapz(dy,trapz(dx,PhiPi(:,:,i).*Di_PhiPi,2))/FL2Norm(i));
    Di_Estimate = I(i)*(FL2Norm(i)*Di_PhiPi + PhiPi(:,:,i).*Di_GPi)/(FL2Norm(i))^2;
    Integrand = Descrepency.*Di_Estimate.*A1;
    GradDescrepency(i+3*SourceNum)  = trapz(dy,trapz(dx,Integrand,2));
    GradPseudo(i+3*SourceNum) = GradDescrepency(i+3*SourceNum) ;
end
end

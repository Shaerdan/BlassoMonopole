function [GradPseudo] = GradPseudoLBFGSB(X,Descrepency,PhiComponent,CosGamma,DistancePQ,Mesh,FL2Norm, Lambda)

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
GradPseudo = zeros(SourceNum,1);
GradSub = GradPseudo;
GradDescrepency  = GradPseudo;

%% Compute the Pseudo gradient of the objective function:
for i = 1:SourceNum
    %% dObjectiveFunction/dIntensity
    integrand_da = (Descrepency.*PhiComponent(:,:,i).*A1);
    GradDescrepency(i)  = trapz(dy,trapz(dx,integrand_da,2));
    GradSub(i) = Lambda*sign(X(i));
    GradPseudo(i) =   GradDescrepency(i)/FL2Norm(i) + GradSub(i);
    
    %% dObjectiveFunction/dRadius 
    Diff_DistPQ_Radius = (Locations(i)-CosGamma(:,:,i))./DistancePQ(:,:,i);
    Diff_Estimate_Radius     = ((-2*I(i)./(DistancePQ(:,:,i).^2)).*Diff_DistPQ_Radius ...
        - (I(i)./(1-Locations(i)*CosGamma(:,:,i)+DistancePQ(:,:,i)))...
        .*(-CosGamma(:,:,i) + Diff_DistPQ_Radius));
    Diff_Estimate_Radius     = I(i)*Diff_Estimate_Radius./FL2Norm(i);
    Diff_FL2Norm_Radius  = - I(i)*((PhiComponent(:,:,i).^2)/(FL2Norm(i))^3).*Diff_Estimate_Radius;
    Diff_Integrand_Radius = Descrepency.*(Diff_Estimate_Radius+Diff_FL2Norm_Radius).*A1;
    GradDescrepency(i+SourceNum)  = trapz(dy,trapz(dx,Diff_Integrand_Radius,2));
    GradPseudo(i+SourceNum) = GradDescrepency(i+SourceNum) ;
    
    %% dObjectiveFunction/dTheta 
    Diff_CosGamma_Theta = -sin(Locations(i+SourceNum))*A2+cos(Locations(i+SourceNum))*A1...
        .*cos(Mesh.PsiQ - Locations(i+2*SourceNum));
    Diff_DistPQ_Theta = -(Locations(i)./DistancePQ(:,:,i)).*(Diff_CosGamma_Theta);
    Diff_Estimate_Theta = ((-2*I(i)./(DistancePQ(:,:,i).^2)).*Diff_DistPQ_Theta ...
        - (I(i)./(1-Locations(i)*CosGamma(:,:,i)+DistancePQ(:,:,i)))...
        .*(-Locations(i)*Diff_CosGamma_Theta + Diff_DistPQ_Theta));
    Diff_Estimate_Theta     = I(i)*Diff_Estimate_Theta./FL2Norm(i);
    Diff_FL2Norm_Theta  = - I(i)*((PhiComponent(:,:,i).^2)/(FL2Norm(i))^3).*Diff_Estimate_Theta;
    Diff_Integrand_Theta = Descrepency.*(Diff_Estimate_Theta+Diff_FL2Norm_Theta).*A1;
    GradDescrepency(i+2*SourceNum) = trapz(dy,trapz(dx,Diff_Integrand_Theta,2));
    GradPseudo(i+2*SourceNum) = GradDescrepency(i+2*SourceNum);
    
    %% dObjectiveFunction/dPsi 
    Diff_CosGamma_Psi =  sin(Locations(i+SourceNum)).*A1...
        .*sin(Mesh.PsiQ - Locations(i+2*SourceNum));
    Diff_DistPQ_Psi = -(Locations(i)./DistancePQ(:,:,i)).*(Diff_CosGamma_Psi);
    Diff_Estimate_Psi   = ((-2*I(i)./(DistancePQ(:,:,i).^2)).*Diff_DistPQ_Psi ...
        -(I(i)./(1-Locations(i)*CosGamma(:,:,i)+DistancePQ(:,:,i)))...
        .*(-Locations(i)*Diff_CosGamma_Psi + Diff_DistPQ_Psi));
    Diff_Estimate_Psi     = I(i)*Diff_Estimate_Psi./FL2Norm(i);
    Diff_FL2Norm_Psi  =- I(i)*((PhiComponent(:,:,i).^2)/(FL2Norm(i))^3).*Diff_Estimate_Psi;
    Diff_Integrand_Psi   = Descrepency.*(Diff_Estimate_Psi+Diff_FL2Norm_Psi).*A1;
    GradDescrepency(i+3*SourceNum) = trapz(dy,trapz(dx,Diff_Integrand_Psi,2));
    GradPseudo(i+3*SourceNum) = GradDescrepency(i+3*SourceNum);
end

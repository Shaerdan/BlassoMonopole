function [GradPseudo] = GradPseudo(X,D,PhiComponent,CosGamma,L_in,Mesh,norm_F, lambda)


SourceNum = length(X)/4;
Intensity = X(1:SourceNum);
Locations = X(SourceNum +1:end);
A1 = sin(Mesh.ThetaQ);
A2 = cos(Mesh.ThetaQ);
dx = Mesh.ThetaQ(1,:)';
dy = Mesh.PsiQ(:,1);
GradPseudo = zeros(1,SourceNum);
GradSub = GradPseudo;
GradDescrepency  = GradPseudo;

for i = 1:SourceNum
    %%%%%%%%%% d/dIntensity %%%%%%%%%%%%%
    integrand_da = (D.*PhiComponent(:,:,i).*A1);
    GradDescrepency(i)  = trapz(dy,trapz(dx,integrand_da,2));
    GradSub(i) = lambda*sign(X(i));
    GradPseudo(i) =   GradDescrepency(i)/norm_F(i) + GradSub(i);
    
    %%%%%%%%%% Delta_r %%%%%%%%%%%%%
    Diff_Lin_Radius = (Locations(i)-CosGamma(:,:,i))./L_in(:,:,i);
    Diff_Estimate_Radius     = ((-2*Intensity(i)./(L_in(:,:,i).^2)).*Diff_Lin_Radius ...
        - (Intensity(i)./(1-Locations(i)*CosGamma(:,:,i)+L_in(:,:,i)))...
        .*(-CosGamma(:,:,i) + Diff_Lin_Radius));
    Diff_Estimate_Radius     = Intensity(i)*Diff_Estimate_Radius./norm_F(i);
    Diff_FL2Norm_Radius  = - Intensity(i)*((PhiComponent(:,:,i).^2)/(norm_F(i))^3).*Diff_Estimate_Radius;
    Diff_Integrand_Radius = D.*(Diff_Estimate_Radius+Diff_FL2Norm_Radius).*A1;
    GradDescrepency(i+SourceNum)  = trapz(dy,trapz(dx,Diff_Integrand_Radius,2));
    GradPseudo(i+SourceNum) = GradDescrepency(i+SourceNum) ;
    
    %%%%%%%%%% Delta_th %%%%%%%%%%%%%
    Diff_CosGamma_Theta = -sin(Locations(i+SourceNum))*A2+cos(Locations(i+SourceNum))*A1...
        .*cos(Mesh.PsiQ - Locations(i+2*SourceNum));
    Diff_Lin_Theta = -(Locations(i)./L_in(:,:,i)).*(Diff_CosGamma_Theta);
    Diff_Estimate_Theta = ((-2*Intensity(i)./(L_in(:,:,i).^2)).*Diff_Lin_Theta ...
        - (Intensity(i)./(1-Locations(i)*CosGamma(:,:,i)+L_in(:,:,i)))...
        .*(-Locations(i)*Diff_CosGamma_Theta + Diff_Lin_Theta));
    Diff_Estimate_Theta     = Intensity(i)*Diff_Estimate_Theta./norm_F(i);
    Diff_FL2Norm_Theta  = - Intensity(i)*((PhiComponent(:,:,i).^2)/(norm_F(i))^3).*Diff_Estimate_Theta;
    Diff_Integrand_Theta = D.*(Diff_Estimate_Theta+Diff_FL2Norm_Theta).*A1;
    GradDescrepency(i+2*SourceNum) = trapz(dy,trapz(dx,Diff_Integrand_Theta,2));
    GradPseudo(i+2*SourceNum) = GradDescrepency(i+2*SourceNum);
    
    %%%%%%%%%% Delta_psi %%%%%%%%%%%%%
    Diff_CosGamma_Psi =  sin(Locations(i+SourceNum)).*A1...
        .*sin(Mesh.PsiQ - Locations(i+2*SourceNum));
    Diff_Lin_Psi = -(Locations(i)./L_in(:,:,i)).*(Diff_CosGamma_Psi);
    Diff_Estimate_Psi   = ((-2*Intensity(i)./(L_in(:,:,i).^2)).*Diff_Lin_Psi ...
        -(Intensity(i)./(1-Locations(i)*CosGamma(:,:,i)+L_in(:,:,i)))...
        .*(-Locations(i)*Diff_CosGamma_Psi + Diff_Lin_Psi));
    Diff_Estimate_Psi     = Intensity(i)*Diff_Estimate_Psi./norm_F(i);
    Diff_FL2Norm_Psi  =- Intensity(i)*((PhiComponent(:,:,i).^2)/(norm_F(i))^3).*Diff_Estimate_Psi;
    Diff_Integrand_Psi   = D.*(Diff_Estimate_Psi+Diff_FL2Norm_Psi).*A1;
    GradDescrepency(i+3*SourceNum) = trapz(dy,trapz(dx,Diff_Integrand_Psi,2));
    GradPseudo(i+3*SourceNum) = GradDescrepency(i+3*SourceNum);
end
GradPseudo = GradPseudo';

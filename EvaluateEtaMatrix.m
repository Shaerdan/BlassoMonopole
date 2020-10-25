function [Eta] = EvaluateEtaMatrix(Lambda,Mesh,P,SwitchNorm)

%   Evaluate Eta over a descrete volumetric mesh:

%   Eta -- ref Equation (45), Page 6 of the report:
%          Eta = Eta(Mesh.S) = integral (
%          F_Normalized(Mesh.Q,Mesh.S)*P(Mesh.Q)dMesh.Q ),
%          dMesh.Q = sin(Mesh.ThetaQ) * dMesh.ThetaQ * dMesh.PsiQ.
%   F  --- The integral kernel, ref Equation (45), Page 6 of the report.
%   FL2Norm --- L2Norm of F;

% dx = Mesh.ThetaQLine;
% dy = Mesh.PsiQLine;
global LogSwitch eps
%%%%%%%%%% testing
%%
ThetaQ = Mesh.ThetaQ;
PsiQ   = Mesh.PsiQ;
SpikeTheta  = Mesh.ThetaSpike(:);
SpikePsi    = Mesh.PsiSpike(:);
SpikeRadius = Mesh.RadiusSpike(:);
N_Spike = length(SpikeTheta);
% tic
% for i = 1 : N_Spike
%     CosGamma = cos(SpikeTheta(i))*cos(ThetaQ(:))+sin(SpikeTheta(i))*...
%         sin(ThetaQ(:)).*cos(SpikePsi(i)- PsiQ(:));
%     DistanceSQ  = sqrt(1 + SpikeRadius(i)^2 - 2* SpikeRadius(i)*CosGamma);
%     Phi(:,i) = (2./DistanceSQ) - log(1 - SpikeRadius(i)*CosGamma + DistanceSQ);
%     if (SwitchNorm == 1)
%         Phi(:,i) = Phi(:,i)/norm(Phi(:,i),2);
%     end
% end
% toc
tic
CosGamma = cos(SpikeTheta').*cos(ThetaQ(:))+sin(SpikeTheta').*...
    sin(ThetaQ(:)).*cos(SpikePsi'- PsiQ(:));
DistanceSQ  = sqrt(1 + SpikeRadius'.^2 - 2* SpikeRadius'.*CosGamma)+ eps;
Phi = (2./DistanceSQ) - LogSwitch*log(1 - SpikeRadius'.*CosGamma + DistanceSQ);

if (SwitchNorm == 1)
    for i=1:N_Spike
        normstore(:,i) = norm(Phi(:,i));
        Phi(:,i) = Phi(:,i)/norm(Phi(:,i));
    end
end
toc
% tic
% [U S V] = svd(Phi);
% toc
% Lambda2 = 1e-10;
% S1 = diag(diag((S.^2+Lambda2^2)./S));
% [n,m] = size(S1);
% S(1:n,1:m) = S1;
% Phi = U*S*V;
PhiS_P = Phi'*P;

% Test = max(P(:))/Lambda;
% PhiS_P2 = Phi3'*P(:);
%     test = Phi - Phi_t
% figure(17)
% surf(reshape(P,size(Mesh.ThetaQ)));
% figure(18)
% plot(PhiS_P2)
% figure(19)
% plot(Phi(:))
% figure(20)
% plot(P(:))

Eta = abs(PhiS_P)/Lambda;


% Eta = abs(Eta)/Lambda;
end


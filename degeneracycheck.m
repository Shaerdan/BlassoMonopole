clc;
clear all;
close all;
global LogSwitch eps;
LogSwitch = 1; eps = 0;
% Linear solver:
% 'FBS' -- linear proximal forward backward splitting;
% 'fminunc' -- matlab solver with a default quasi-newton algorithm
% and bfgs approximation of the Hessian.
% Set the nonlinear solver:
% 'lbfgsc' -- limited memory bfgs solver.
% 'fminunc' -- matlab solver with a default quasi-newton algorithm
%  and bfgs approximation of the Hessian.
%% Setup Parameters:
%%% Regularization
Lambda = 0.2*1e-16; GlobalIteration = 10; IntensityMax = 20;
eps_Eta = 1e-1;
% Lambda = .2*1e-5;
% Lambda = .5*1e-9;
%%% Mesh Control
MeshPointSSize.Radius = 10; MeshPointSSize.Theta = 20; MeshPointSSize.Psi = 40;
MeshPointQSize.Theta = 50; MeshPointQSize.Psi = 50; Cap.RadiusUpper = 0.9999; Cap.RadiusLow = 0.0001;
Cap.Theta = pi; Cap.Psi = 2*pi;
%%  Swtiches
Switch.Normalize = 1;            % Normalization switch: 0 -- disable  1 -- enable.
Switch.Max = 'grid';             % Maximization method: gradient -- fmincon,
% grid -- discrete grid,  fixed -- fixed location
% gridgrad -- grid + fmincon
Switch.Display = 1;              % Display figures in each loop 1 or not 0
Switch.LinearSolver = 'admm';             % choose 'fbs' or 'admm' or 'none';
Switch.NonLinearSolver = 'lbfgsc';        % choose 'lbfgsc' or 'fmincon';
Swtich.Conditioning = 0;
%% Exact Souce Parameters:
RadiusExact = [0.3 0.3]; SourceNum = 2; % Exact source depths and source number setup.
IntensityExact = [1 -1]; % Exact source intensity setup.
ThetaExact = [ 1.5 1.501 ]; PsiExact = [ 1.5 1.501 ]; % Exact source angle setup.
NoiseLevel = 0;
LocationExact = zeros(3*SourceNum,1);   % Initialize real source location assembly vector.
LocationExact (1:SourceNum) = RadiusExact; % Store radius component to the assembly vector.
LocationExact (SourceNum+1:2*SourceNum) = ThetaExact; % Store theta component ....
LocationExact (2*SourceNum+1:end) = PsiExact; % Store psi component ...
Mesh = GenerateMesh(MeshPointSSize,MeshPointQSize,Cap);

ThetaQ = Mesh.ThetaQ;
PsiQ   = Mesh.PsiQ;
SpikeTheta  = Mesh.ThetaSpike(:);
SpikePsi    = Mesh.PsiSpike(:);
SpikeRadius = Mesh.RadiusSpike(:);
SpikeTheta0  = ThetaExact';
SpikePsi0    = PsiExact';
SpikeRadius0 = RadiusExact';
N_Spike = length(SpikeTheta);

tic
CosGamma = cos(SpikeTheta').*cos(ThetaQ(:))+sin(SpikeTheta').*...
    sin(ThetaQ(:)).*cos(SpikePsi'- PsiQ(:));
DistanceSQ  = sqrt(1 + SpikeRadius'.^2 - 2* SpikeRadius'.*CosGamma)+ eps;
Phi = (2./DistanceSQ) - LogSwitch*log(1 - SpikeRadius'.*CosGamma + DistanceSQ);

CosGamma0 = cos(SpikeTheta0').*cos(ThetaQ(:))+sin(SpikeTheta0').*...
    sin(ThetaQ(:)).*cos(SpikePsi0'- PsiQ(:));
DistanceSQ0  = sqrt(1 + SpikeRadius0'.^2 - 2* SpikeRadius0'.*CosGamma0)+ eps;
Phi0 = (2./DistanceSQ0) - LogSwitch*log(1 - SpikeRadius0'.*CosGamma0 + DistanceSQ0);
if (Switch.Normalize == 1)
    for i=1:N_Spike
        normstore(:,i) = norm(Phi(:,i)); 
        Phi(:,i) = Phi(:,i)/norm(Phi(:,i));
    end
    for i=1:2
        normstore0(:,i) = norm(Phi0(:,i)); 
        Phi0(:,i) = Phi0(:,i)/norm(Phi0(:,i)); 
    end
    
end
toc
X = [SpikeRadius0',SpikeTheta0',SpikePsi0'];
tic
dPhi0 = Compute_dPhi(X,Mesh,Switch.Normalize);
toc
GammaX = [Phi0,dPhi0];
tic
GammaPlus = GammaX'*GammaX;
toc
tic
CondGammaPlus = cond(GammaPlus);
RankGammaPlus = rank(CondGammaPlus);
toc
Gamma2 = GammaX*GammaX';
CondGamma2 = cond(Gamma2);
GammaPlusInv = inv(GammaPlus);


eta_v = Phi'*GammaX*(inv(GammaPlus)*[-1 1 0 0 0 0 0 0]');
max(eta_v)

% Test = Phi(:,1)-Phi(:,2);
% Lambda2 = 1e1;
% [U S V] = svd(Phi);
%     S1 = diag(diag((S.^2+Lambda2^2)./S));
%     [n,m] = size(S1);
%     S(1:n,1:m) = S1;
%     Phi2 = U*S*V;
%     CondPhi = cond(Phi'*Phi); CondPhi2 = cond(Phi2'*Phi2);
% Phi = (A + Lambda2*eye(size(A)));
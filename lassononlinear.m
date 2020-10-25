clc;
clear all;
close all;
%% Setup Parameters:
%%% Regularization
Lambda = 1e-8;

%%% Mesh Control
MeshPointSSize.Radius = 10; MeshPointSSize.Theta = 10; MeshPointSSize.Psi = 10;
MeshPointQSize.Theta = 100; MeshPointQSize.Psi = 100; Cap.Radius = 0.9;
Cap.Theta = pi; Cap.Psi = 2*pi;

%% Souce Control:
RadiusReal = [0.8 0.8]; SourceNum = 2; % Real source depths and source number setup.
IntensityReal = [1 -1]; % Real source intensity setup.
ThetaReal = [.2 .25 ]; PsiReal = [.2 .25 ]; % Real source angle setup.
NoiseLevel = 1e-2;
LocationReal = zeros(3*SourceNum,1);   % Initialize real source location assembly vector.
LocationReal (1:SourceNum) = RadiusReal; % Store radius component to the assembly vector.
LocationReal (SourceNum+1:2*SourceNum) = ThetaReal; % Store theta component ....
LocationReal (2*SourceNum+1:end) = PsiReal; % Store psi component ...

%% Create the mesh for point Q in R^3 and point S in R^2 (Q,S as defined
% in the report).
% Mesh:
% Mesh.RadiusQ -- Radius component of Q mesh.
% Mesh.ThetaQ -- Theta component of Q mesh.
% Mesh.PsiQ -- Psi component of Q mesh.
% Mesh.ThetaS -- Theta component of S mesh.
% Mesh.PsiS -- Psi component of S mesh.
Mesh = GenerateMesh(MeshPointSSize,MeshPointQSize,Cap);

%% Create measurement data.

% Data: Data.Measurement, Data.Noise
Data = GenerateData(IntensityReal,LocationReal,Mesh,NoiseLevel);

%% Choose method:
method = 'lbfgsc';       % 'fmincon' or 'lbfgsc'

%% fminunc solver control:
Opts = optimoptions(@fmincon,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',10000,'MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-30,'CheckGradients',true,...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',false);


%% Initializing %%%%%%

% InitialSolution = 10*randn(4*SourceNum,1);
InitialSolution = zeros(4*SourceNum,1);
IntensityMax = 50;
LowerBoundNonLinear = [-IntensityMax*ones(1,SourceNum), ...
    repmat(0.001,1,SourceNum),zeros(1,SourceNum),...
    zeros(1,SourceNum)];
UpperBoundNonLinear = [IntensityMax*ones(1,SourceNum),...
    repmat(Cap.Radius,1,SourceNum),repmat(pi,1,SourceNum),...
    repmat(2*pi,1,SourceNum)];
% LowerBoundNonLinear = []; UpperBoundNonLinear = [];
switch method
    case 'fmincon'
        problem = createOptimProblem('fmincon','x0',InitialSolution,...
            'ub',UpperBoundNonLinear,'lb',LowerBoundNonLinear, ...
            'objective',@(X) ObjectiveFuncNonLinearFminunc(X,Data.Measurement,...
            Lambda, Mesh),...
            'options', Opts);
        [SolutionArgmin,fval,exitflag,output,lambda,grad,hessian] = fmincon(problem);
%         fprintf('flag of fmincon %d\n',flag);
        Solution.Intensity = SolutionArgmin(1:SourceNum)';
        Solution.Location = SolutionArgmin(SourceNum+1:end)';
    case 'lbfgsc'
        OptLBFGSB = struct('x0',InitialSolution,'m',100,'factr',1e0,...
            'pgtol',1e-16,'maxIts',5000,'maxTotalIts',50000,'printEvery',1);
        [SolutionArgmin,f_val2,info2] = lbfgsb( @(X) ...
            ObjectiveFuncNonLinearLBFGSB(X,Data.Measurement,Lambda, Mesh),...
            LowerBoundNonLinear', UpperBoundNonLinear', OptLBFGSB );
        Solution.Intensity = SolutionArgmin(1:SourceNum)';
        %             if (ii>4)
        %                 Solution.Intensity(abs(Solution.Intensity)<= Cropping) = 0;
        %             end
        Solution.Location = SolutionArgmin(SourceNum+1:end)';
end

IntensitySoln = Solution.Intensity;
figure(4)
bar(IntensitySoln);
RadiusSoln = Solution.Location(1:SourceNum);
ThetaSoln = Solution.Location(SourceNum+1:2*SourceNum);
PsiSoln = Solution.Location(2*SourceNum+1:3*SourceNum);
figure(6)
CartesianXSoln = RadiusSoln.*sin(ThetaSoln).*cos(PsiSoln);
CartesianYSoln = RadiusSoln.*sin(ThetaSoln).*sin(PsiSoln);
CartesianZSoln = RadiusSoln.*cos(ThetaSoln);
CartesianXReal = RadiusReal.*sin(ThetaReal).*cos(PsiReal);
CartesianYReal = RadiusReal.*sin(ThetaReal).*sin(PsiReal);
CartesianZReal = RadiusReal.*cos(ThetaReal);
[XSphere,YSphere,ZSphere] = sphere;
h = surf(XSphere, YSphere, ZSphere,'FaceColor','none');
set(h, 'FaceAlpha', 0.2);
hold on;
for i = 1:SourceNum
    scatter3(CartesianXReal(i),CartesianYReal(i),CartesianZReal(i),100,'r','linewidth',1.5);
end
for i = 1:SourceNum
    scatter3(CartesianXSoln(i),CartesianYSoln(i),CartesianZSoln(i),100,'k','*','linewidth',1.5);
    hold on;
end
hold off;



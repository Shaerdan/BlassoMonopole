clc;
clear all;
close all;

%% Setup Parameters:
%%% Regularization
Lambda = 1e-5;

%%% Mesh Control
meshsize_r = 30; meshsize_theta = 30; meshsize_psi = 30;

%% Souce Control:
RadiusReal = 0.5; SourceNum = 2; % Real source depths and source number setup.
IntensityReal = [1 -1]; % Real source intensity setup.
ThetaReal = [0.4, 0.5]; PsiReal = [1.2, 1.3]; % Real source angle setup.
NoiseLevel = 0;
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
Mesh = GenerateMesh(meshsize_r,meshsize_theta,meshsize_psi);

%% Create measurement data.

% Data: Data.Measurement, Data.Noise
Data = GenerateData(IntensityReal,LocationReal,Mesh,NoiseLevel);

%% Choose solvers:
% Linear solver:
% 'FBS' -- linear proximal forward backward splitting;
% 'fminunc' -- matlab solver with a default quasi-newton algorithm
% and bfgs approximation of the Hessian.
LinearSolver = 'FBS';               % choose 'FBS' or 'fminunc';
% Set the nonlinear solver:
% 'lbfgsc' -- limited memory bfgs solver.
% 'fminunc' -- matlab solver with a default quasi-newton algorithm
%  and bfgs approximation of the Hessian.
NonLinearSolver = 'fminunc';        % choose 'lbfgsc' or 'fminunc';

%% FBS solver control:
FBS.Tau = 1e-1;
FBS.GradTol = 1e-16; FBS.StepTol = 1e-16;
FBS.Iteration = 5e5;
FBS.DisplayFrequency = 100;

%% fminunc solver control:
Opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',10000,'MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-10,'StepTolerance',1e-16);


%% Initializing the Sliding Frank Wolfe iteration %%%%%%
GlobalIteration = 1000;

Solution.Intensity= ([]); Solution.Location = ([]);
% Solution.Radius .Theta .Psi below is used to update the Location assembly vector.
Solution.Radius = ([]); Solution.Theta = ([]); Solution.Psi = ([]);
% Initialize EtaMax.
EtaMax = zeros(1,GlobalIteration);


%% Sliding Frank Wolfe Algorithm Interation:
for ii = 1:GlobalIteration
    
    F_L2Norm = L2NormF(Solution.Location,Mesh);
    
    %%  Solve (r,theta,psi) = Argmax_(r',theta',psi') {Eta(r',theta',psi')}
    %          Eta = Eta(Mesh.S) = integral (
    %          F_Normalized(Mesh.Q,Mesh.S)*P(Mesh.Q)dMesh.Q ),
    %          dMesh.Q = sin(Mesh.ThetaQ) * dMesh.ThetaQ * dMesh.PsiQ.
    P   = EvaluateP(Solution.Location,Solution.Intensity,Data.Measurement,Mesh,F_L2Norm);
    %%% MARKED FOR DEBUGGING %%%
    tic
    Eta = EvaluateEta(Lambda,Mesh,P);
    toc
    [EtaMax(ii),Index] = max(Eta(:));
    [IndexRadius,IndexTheta,IndexPsi] = ind2sub(size(Eta),Index);
    %%% Debug line, output EtaMax:
    fprintf('Eta= %f', EtaMax(ii));
    %%% MARKED FOR DEBUGGING %%%
    
    % obtain the new source location from the index of EtaMax.
    RadiusUpdate = Mesh.RadiusS(IndexRadius);
    ThetaUpdate = Mesh.ThetaS(IndexTheta);
    PsiUpdate = Mesh.PsiS(IndexPsi);
    
    figure(2)
    bar([RadiusUpdate, ThetaUpdate, PsiUpdate]);
    
    %%% Termination condition for the global interation
    if (EtaMax(ii) <=1)
        break
    end
    
    % Update the source location assembly vector:
    Solution.Radius = [Solution.Radius,RadiusUpdate];
    Solution.Theta = [Solution.Theta, ThetaUpdate];
    Solution.Psi = [Solution.Psi, PsiUpdate];
    Solution.Location = [Solution.Radius, Solution.Theta, Solution.Psi];
    % End of (r,theta,psi) = Argmax_(r',theta',psi'){Eta(r',theta',psi')} Procedure
    
    %% Update the source intensity by solving the linear lasso problem
    
    % Initialize the linear minimization problem:
    IntensityUpdate = 0;
    InitialIntensity = [Solution.Intensity, IntensityUpdate];
    PhiComponent = ComputePotentialComponent(Solution.Location,Mesh);
    F_L2Norm = L2NormF(Solution.Location,Mesh);
    % Choose solvers:
    switch LinearSolver
        case 'FBS'
            [Solution.Intensity] = FBSSolver(InitialIntensity,Data.Measurement,...
                FBS.Iteration,PhiComponent,Mesh,...
                FBS.GradTol,FBS.StepTol,FBS.Tau,Lambda,...
                FBS.DisplayFrequency,F_L2Norm);
        case 'fminunc'
            SourceNumUpdate = ii;
            IntensityMax = 5;
            % Lower bound = -IntensityMax* RealIntensity:
            LowerBoundLinear = -IntensityMax*ones(1,SourceNumUpdate)';
            % Upper bound =  IntensityMax* RealIntensity
            UpperBoundLinear = IntensityMax*ones(1,SourceNumUpdate)';
            problem = createOptimProblem('fmincon','x0',InitialIntensity,...
                'ub',UpperBoundLinear,'lb',LowerBoundLinear, ...
                'objective',@(X) ObjectiveFuncLinearLasso(X,...
                Data.Measurement,PhiComponent,Mesh,F_L2Norm,Lambda),...
                'options',Opts);
            [Solution.Intensity] = fmincon(problem);
    end
    
    % Debug line, output intensities:
    figure(22)
    bar(Solution.Intensity);
    % End of source intensity update
    
    %% Update source intensity and locations by solving anonlinear Lasso problem:
    % Initialize the nonlinear minimization problem:
    InitialSolution = [Solution.Intensity,Solution.Location];
    SourceNumUpdate = length(Solution.Intensity);
    % The box constraint is set to be: -IntensityMax < Intensity <
    % IntensityMax, 0 < Radius < 0.999, 0 < Theta < pi, 0 < Psi < 2*pi:
    IntensityMax = 5;
    LowerBoundNonLinear = [-IntensityMax*ones(1,SourceNumUpdate), ...
        repmat(0.001,1,SourceNumUpdate),zeros(1,SourceNumUpdate),...
        zeros(1,SourceNumUpdate)];
    UpperBoundNonLinear = [IntensityMax*ones(1,SourceNumUpdate),...
        repmat(0.999,1,SourceNumUpdate),repmat(pi,1,SourceNumUpdate),...
        repmat(2*pi,1,SourceNumUpdate)];
    
    % Choose nonlinear solver:
    switch NonLinearSolver
        case 'lbfgsc'
            options1 = struct('x0',InitialSolution,'m',20,'factr',1e0,...
                'pgtol',1e-30,'maxIts',10000,'maxTotalIts',500000,'printEvery',1);
            IntensityMax = 5;
            [SolutionArgmin,f_val2,info2] = lbfgsb( @(X) ...
                ObjectiveFuncNonLinearLasso(X,Data,lambda, Mesh, F_L2Norm),...
                lb1, ub1, options1 );
        case 'fminunc'
            problem = createOptimProblem('fmincon','x0',InitialSolution,...
                'ub',UpperBoundNonLinear,'lb',LowerBoundNonLinear, ...
                'objective',@(X) ObjectiveFuncNonLinearLasso(X,Data.Measurement,...
                Lambda, Mesh,F_L2Norm),...
                'options', Opts);
            [SolutionArgmin,fval2,flag,stepcount] = fmincon(problem);
            fprintf('flag of fmincon %d\n',flag);
    end
    Solution.Intensity = SolutionArgmin(1:SourceNumUpdate);
    Solution.Location = SolutionArgmin(SourceNumUpdate+1:end);
    % Debug line, output the source intensities and locations:
    figure(23)
    bar(Solution.Intensity);
    figure(24)
    bar(Solution.Location);
    % End of source locations and intensities update.
    
    %% Display The Solution:
    IntensitySoln = Solution.Intensity;
    RadiusSoln = Solution.Location(1:SourceNumUpdate);
    ThetaSoln = Solution.Location(SourceNumUpdate+1:2*SourceNumUpdate);
    PsiSoln = Solution.Location(2*SourceNumUpdate+1:3*SourceNumUpdate);
    figure(11)
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
    for i = 1:SourceNumUpdate
        scatter3(CartesianXSoln(i),CartesianYSoln(i),CartesianZSoln(i),100,'k','*','linewidth',1.5);
        hold on;
    end
    
end






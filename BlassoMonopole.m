clc;
clear all;
close all;
global LogSwitch eps;
LogSwitch = 1; eps = 1e-3;
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
MeshPointQSize.Theta = 50; MeshPointQSize.Psi = 50; Cap.RadiusUpper = 0.9999;Cap.RadiusLow = 0.0001;
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
Data = GenerateData(IntensityExact,LocationExact,Mesh,NoiseLevel,Switch.Normalize);
%% Choose solvers:


%% FBS solver control:
FBS.GradTol = 1e-16; FBS.StepTol = 1e-16;
FBS.Iteration = 5e5;
FBS.DisplayFrequency = 100000;
FBS.Mu = 0.4;

%% ADMM solver control:
ADMM.Tol = 1e-16; ADMM.Steps = 1000; ADMM.Rho = 0.5;

%% fminunc solver control:
Opts = optimoptions(@fmincon,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',10000,'MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-16);

% Opts = optimoptions(@fmincon,...
%     'MaxIterations',10000,'MaxFunctionEvaluations',20000,...
%     'OptimalityTolerance',1e-16,'StepTolerance',1e-16);

Opts1 = optimoptions(@fmincon,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',10000,'MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-16,'CheckGradients',false,...
    'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false);

%% Initializing the Sliding Frank Wolfe iteration %%%%%%
initiation = 'cold';         InitialSize = 0;

switch initiation
    case 'cold'
        Solution.Intensity= ([]); Solution.Location = ([]);
        % Solution.Radius .Theta .Psi below is used to update the Location assembly vector.
        Solution.Radius = ([]); Solution.Theta = ([]); Solution.Psi = ([]);
        InitialSize = 0;
    case 'warm'
        Solution.Intensity= rand(1,InitialSize);
        Solution.Radius = 0.99*rand(1,InitialSize); Solution.Theta = pi*rand(1,InitialSize);
        Solution.Psi = 2*pi*rand(1,InitialSize);
        Solution.Location = [Solution.Radius,Solution.Theta,Solution.Psi];
end
EtaMax = zeros(1,GlobalIteration);
IntensityLabels = cell(GlobalIteration+1,1);
for i = 1:GlobalIteration+1
    IntensityLabels{i} = sprintf('I%d',i);
end

% Set cropping threshold:
Cropping = 1e-1;
UpdatesRecording = [];
%% Sliding Frank Wolfe Algorithm Interation:
for ii = InitialSize+1:GlobalIteration+InitialSize+1
    %     if (Switch.Normalize == 1)
    %         FL2Norm = L2NormF(Solution.Location,Mesh);
    %     else
    %         FL2Norm = 1;
    %     end
    
    %%  Solve (r,theta,psi) = Argmax_(r',theta',psi') {Eta(r',theta',psi')}
    
    P   = EvaluateP(Solution.Location,Solution.Intensity,Data.Measurement,Mesh,Switch.Normalize);
    %%% MARKED FOR DEBUGGING %%%
    switch Switch.Max
        case 'grid'
            Mesh = GenerateMesh(MeshPointSSize,MeshPointQSize,Cap);
            Eta = EvaluateEtaMatrix(Lambda,Mesh,P,Switch.Normalize);
            [EtaMax(ii),Index] = max(Eta);
            %     %%% Debug line, output EtaMax:
            fprintf('Eta= %f', EtaMax(ii));
            
            %     pause(1);
            SpikeTheta  = Mesh.ThetaSpike(:);
            SpikePsi    = Mesh.PsiSpike(:);
            SpikeRadius = Mesh.RadiusSpike(:);
            RadiusUpdate = SpikeRadius(Index);
            ThetaUpdate = SpikeTheta(Index);
            PsiUpdate = SpikePsi(Index);
            UpdatesRecording = [UpdatesRecording, RadiusUpdate,ThetaUpdate,PsiUpdate];
            if (Switch.Display == 1)
                figure(2)
                Yaxs = [RadiusUpdate ThetaUpdate PsiUpdate];
                Xaxs = categorical({'Radius','Theta','Psi'});
                Xaxs = reordercats(Xaxs,{'Radius','Theta','Psi'});
                bar(Xaxs,Yaxs);
                legend('New Spike Location from Argmax(\eta_k)')
            end
            if (Switch.Display == 1)
                figure(20)
                plot(log10(EtaMax));
                xlabel('Iteration, k');
                ylabel('log_{10}(|\eta_k|_{max})');
                legend('Convergence history of max|\eta_k|')
            end
            if (abs(EtaMax(ii)) <=1+eps_Eta)
                break
            end
        case 'gradient'
            %%  Gradient ascent for argmax (-0.5*|eta|^2)
            LowerBoundRefine = [0, 0, 0]';
            UpperBoundRefine = [Cap.RadiusUpper, Cap.Theta, Cap.Psi]';
            InitialSolution = zeros(3,1);
            problem = createOptimProblem('fmincon','x0',InitialSolution,...
                'ub',UpperBoundRefine,'lb',LowerBoundRefine, ...
                'objective',@(X) ObjEtaMatrix(X,Lambda,Mesh,Solution,Data,FL2Norm,Switch.Normalize),...
                'options', Opts);
            [SolutionLoc,EtaMin] = fmincon(problem);
            RadiusUpdate = SolutionLoc(1);
            ThetaUpdate = SolutionLoc(2);
            PsiUpdate = SolutionLoc(3);
            %             figure(2)
            %             bar([RadiusUpdate, ThetaUpdate, PsiUpdate]);
            EtaMax(ii) = 2*sqrt(abs(EtaMin));
            fprintf('Eta= %f', EtaMax(ii));
            
            if (abs(EtaMax(ii)) <=1)
                break
            end
        case 'gridgrad'
            %%  Hybrid Method
            % Coarse maximization:
            Mesh = GenerateMesh(MeshPointSSize,MeshPointQSize,Cap);
            tic
            Eta = EvaluateEtaMatrix(Lambda,Mesh,Solution,Data,Switch.Normalize);
            toc
            [EtaMax(ii),Index] = max(Eta);
            
            %%% Debug line, output EtaMax:
            fprintf('Eta= %f', EtaMax(ii));
            
            SpikeTheta  = Mesh.ThetaSpike(:);
            SpikePsi    = Mesh.PsiSpike(:);
            SpikeRadius = Mesh.RadiusSpike(:);
            RadiusUpdate = SpikeRadius(Index);
            ThetaUpdate = SpikeTheta(Index);
            PsiUpdate = SpikePsi(Index);
            
            % Refine Steps: gradient ascent for argmin[-(1/2)*Eta^2]:
            LowerBoundRefine = [0, 0, 0]';
            UpperBoundRefine = [Cap.RadiusUpper, Cap.Theta, Cap.Psi]';
            InitialSolution = [RadiusUpdate,ThetaUpdate,PsiUpdate]';
            problem = createOptimProblem('fmincon','x0',InitialSolution,...
                'ub',UpperBoundRefine,'lb',LowerBoundRefine, ...
                'objective',@(X) ObjEtaMatrix(X,Lambda, Mesh,Switch.Normalize),...
                'options', Opts);
            [SolutionLoc,EtaMin] = fmincon(problem);
            RadiusUpdate = SolutionLoc(1);
            ThetaUpdate = SolutionLoc(2);
            PsiUpdate = SolutionLoc(3);
            UpdatesRecording = [UpdatesRecording, RadiusUpdate,ThetaUpdate,PsiUpdate];
            
            EtaMax(ii) = 2*sqrt(abs(EtaMin));
            fprintf('Eta= %f', EtaMax(ii));
            if (Switch.Display == 1)
                figure(20)
                plot(log10(EtaMax));
                xlabel('Iteration, k');
                ylabel('log_{10}(|\eta_k|_{max})');
                legend('Convergence history of max|\eta_k|')
            end
            if (abs(EtaMax(ii)) <=1)
                break
            end
        case 'fixed'
            [RadiusUpdate,ThetaUpdate,PsiUpdate] = deal(0.7,pi*(rand),2*pi*(rand));
            EtaMax = 0;
    end
    
    
    
    % Update the source location assembly vector:
    Solution.Radius = [Solution.Radius, RadiusUpdate];
    Solution.Theta = [Solution.Theta, ThetaUpdate];
    Solution.Psi = [Solution.Psi, PsiUpdate];
    Solution.Location = [Solution.Radius, Solution.Theta, Solution.Psi];
    % End of (r,theta,psi) = Argmax_(r',theta',psi'){Eta(r',theta',psi')} Procedure
    
    %% Update the source intensity by solving the linear lasso problem
    
    % Initialize the linear minimization problem:
    IntensityUpdate = 0;
    InitialIntensity = [Solution.Intensity, IntensityUpdate];
    PhiComponent = ComputePotentialComponent(Solution.Location,Mesh,Switch.Normalize);
    %     if (Switch.Normalize == 1)
    %         FL2Norm = L2NormF(Solution.Location,Mesh);
    %     else
    %         FL2Norm = 1;
    %     end
    %% S2 of Algorithm1, the linear lasso problem:
    switch Switch.LinearSolver
        case 'fbs'
            [M,b]= ComputeHessian(Data.Measurement,PhiComponent,FL2Norm,Mesh,...
                Switch.Normalize);
            LipschitzConst = norm(M,2);
            FBS.Tau = 1/LipschitzConst;
            Solution.Intensity = FBSSolver(InitialIntensity',M,b,Data.Measurement,...
                FBS,PhiComponent,Mesh,Lambda,FL2Norm');
            Solution.Intensity = Solution.Intensity';
        case 'admm'
            [M,b]= ComputeHessian(Data.Measurement,PhiComponent,Swtich.Conditioning);
            Solution.Intensity = admm(InitialIntensity',M,b,Lambda,ADMM);
            Solution.Intensity = Solution.Intensity';
        case 'none'
            Solution.Intensity = InitialIntensity;
    end
    
    % Debug line, output intensities:
    if (Switch.Display == 1)
        figure(3)
        Ilabel = categorical(IntensityLabels(1:ii));
        Ilabel = reordercats(Ilabel,IntensityLabels(1:ii));
        bar(Ilabel,Solution.Intensity);
        legend('Intensity Solution from Linear Lasso')
    end
    
    %% Update source intensity and locations by solving anonlinear Lasso problem:
    % Initialize the nonlinear minimization problem:
    InitialSolution = [Solution.Intensity,Solution.Location];
    SourceNumUpdate = length(Solution.Intensity);
    % The box constraint is set to be: -IntensityMax < Intensity <
    % IntensityMax, 0 < Radius < 0.999, 0 < Theta < pi, 0 < Psi < 2*pi:
    LowerBoundNonLinear = [-IntensityMax*ones(1,SourceNumUpdate), ...
        zeros(1,SourceNumUpdate),zeros(1,SourceNumUpdate),...
        zeros(1,SourceNumUpdate)];
    UpperBoundNonLinear = [IntensityMax*ones(1,SourceNumUpdate),...
        repmat(Cap.RadiusUpper,1,SourceNumUpdate),repmat(pi,1,SourceNumUpdate),...
        repmat(2*pi,1,SourceNumUpdate)];
    
    %% S3 of Algorithm 1, the nonlinear lasso problem:
    switch Switch.NonLinearSolver
        case 'lbfgsc'
            OptLBFGSB = struct('x0',InitialSolution','m',100,'factr',1e-16,...
                'pgtol',1e-32,'maxIts',100000,'maxTotalIts',50000000,'printEvery',1);
            [SolutionArgmin,f_val2,info2] = lbfgsb( @(X) ...
                ObjectiveFuncNonLinearLBFGSB(X,Data.Measurement,Lambda, Mesh,Switch.Normalize,Swtich.Conditioning),...
                LowerBoundNonLinear', UpperBoundNonLinear', OptLBFGSB );
            Solution.Intensity = SolutionArgmin(1:SourceNumUpdate)';
            Solution.Location = SolutionArgmin(SourceNumUpdate+1:end)';
        case 'fmincon'
            problem = createOptimProblem('fmincon','x0',InitialSolution',...
                'ub',UpperBoundNonLinear,'lb',LowerBoundNonLinear, ...
                'objective',@(X) ObjectiveFuncNonLinearFminunc(X,Data.Measurement,...
                Lambda, Mesh,Switch.Normalize),...
                'options', Opts1);
            [SolutionArgmin,fval2,flag,stepcount] = fmincon(problem);
            fprintf('flag of fmincon %d\n',flag);
            Solution.Intensity = SolutionArgmin(1:SourceNumUpdate)';
            Solution.Location = SolutionArgmin(SourceNumUpdate+1:end)';
    end
    
    % Update Solution Location Components:
    Solution.Radius = Solution.Location(1:ii);
    Solution.Theta  = Solution.Location(ii+1:2*ii);
    Solution.Psi    = Solution.Location(2*ii+1:3*ii);
    % Debug line, output the source intensities and locations:
    if (Switch.Display == 1)
        figure(4)
        Ilabel = categorical(IntensityLabels(1:ii));
        Ilabel = reordercats(Ilabel,IntensityLabels(1:ii));
        bar(Ilabel,Solution.Intensity);
        legend('Intensity Solution from Nonlinear Lasso')
    end
    
    %% Display The Solution:
    IntensitySoln = Solution.Intensity;
    RadiusSoln = Solution.Location(1:SourceNumUpdate);
    ThetaSoln = Solution.Location(SourceNumUpdate+1:2*SourceNumUpdate);
    PsiSoln = Solution.Location(2*SourceNumUpdate+1:3*SourceNumUpdate);
    if (Switch.Display == 1)
        
        figure(6)
        CartesianXSoln = RadiusSoln.*sin(ThetaSoln).*cos(PsiSoln);
        CartesianYSoln = RadiusSoln.*sin(ThetaSoln).*sin(PsiSoln);
        CartesianZSoln = RadiusSoln.*cos(ThetaSoln);
        CartesianXExact = RadiusExact.*sin(ThetaExact).*cos(PsiExact);
        CartesianYExact = RadiusExact.*sin(ThetaExact).*sin(PsiExact);
        CartesianZExact = RadiusExact.*cos(ThetaExact);
        [XSphere,YSphere,ZSphere] = sphere;
        h = surf(XSphere, YSphere, ZSphere,'FaceColor','none');
        set(h, 'FaceAlpha', 0.2);
        hold on;
        for i = 1:SourceNum
            scatter3(CartesianXExact(i),CartesianYExact(i),CartesianZExact(i),100,'r','linewidth',1.5);
        end
        for i = 1:SourceNumUpdate
            scatter3(CartesianXSoln(i),CartesianYSoln(i),CartesianZSoln(i),100,'k','*','linewidth',1.5);
            hold on;
        end
        hold off;
    end
    
end






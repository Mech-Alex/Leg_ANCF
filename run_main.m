% Multibody dynamics model of human lower limb with ANCF Achilles tendon.
% Based on the code by Erik Makkonen 2021 and Leonid Obrezkov 2022
% Aleksandr Nemov

 clc; close all; 
 clear variables;

addpath('MBS_functions','Tendon_functions')
addpath("Postprocessing");

format long
% inpuit is controlled in eomLeg.m file

% Parameters
%g = 9.81;   % m/s2
g = 0;
% Parameters are from OpenSim gait 2354 model

% Body 1 - femur 
m1 = 9.3014e4;    % kg
L11 = 0.17;     % m, distance to hip joint
%L11 = 0.207;     % m, distance to hip joint
L12 = 0.2258;   % m, distance to knee joint 
%L12 = 0.222;   % m, distance to knee joint 
% Inertia tensors
Ix1_ = 0.1412;
Iy1_ = 0.0351;
Iz1_ = 0.1339;

% Body 2 - tibia and fibula
m2 = 3.7075e4;    % kg
L21 = 0.1867;   % m, distance to knee joint
%L21 = 0.160;   % m, distance to knee joint
L22 = 0.2433;   % m, distance to ankle joint
%L22 = 0.200;   % m, distance to ankle joint
% Inertia tensors
Ix2_ = 0.0511;
Iy2_ = 0.0051;
Iz2_ = 0.0504;

% Body 3 - foot
% gait 2354 has 3 bodies, talus, calcaneum and toes
% that are combined into one body in this model
m3 = 1.5666;    % kg
%L31 = 0.0636;   % m, length distance to ankle joint
%L31 = 0.0533897;   % m, length distance to ankle joint
L31 = 0.091;
%H3 = 0.0148;    % m, height distance to ankle joint
H3 = 0.02;    % m, height distance to ankle joint
% Inertia tensors
Ix3_ = 0.0052;
Iy3_ = 0.0025;
Iz3_ = 0.0051;

% muscle parameters
muscle_sol = muscle_param_sol();
muscle_ta = muscle_param_ta();

% SOL attachment points
% location in local coordinate system
u_sol1 = [0; -0.1126; 0.0038];  % heel attachment, body 3
u_sol2 = [0; -0.002; 0.0337];   % tibia attachment, body 2

% TA attachment points
u_ta1 = [0; 0.018; 0.0247];     % tibia attachment, body 2
u_ta2 = [0; 0.033; -0.2083];    % tibia lower point, body 2
u_ta3 = [0; 0.0044; -0.0092];              % foot, body 3

parameters = [g, m1, L11, L12, Ix1_, Iy1_, Iz1_, m2, L21, L22, ...
    Ix2_, Iy2_, Iz2_, m3, L31, H3, Ix3_, Iy3_, Iz3_, u_sol1.', u_sol2.',...
    u_ta1.', u_ta2.', u_ta3.'];


%%%
% generate symbolic MBS equations
if (exist('body')==0)
body = generateLeg();
end

%%%

% Initial conditions
y0 = zeros(77,1);               % 42 coordinates 3+5+5+5 constraints, 
                                % 8 muscle parameters
                                % 8 for wanted results
y0(3) = -L11;                   % body1 z coordinate at start
y0(4) = 1;                      % e01
y0(10) = y0(3)-L12-L21;         % body2 z loc
y0(11) = 1;                     % e02
y0(16) = L31;                   % body3 y loc
y0(17) = y0(10)-L22-H3;         % body3 z loc
y0(18) = 1;                     % e03

% muscles
% sol initial mtc length
% distance between 2 muscle attachment points 
y0(61) = l_mtc_2points(y0(15:21), y0(8:14), u_sol1, u_sol2);   

% CE length
q_sol = 0.001;  % initial activation level [0 < q =< 1]
% 0.001 is minimal activation, no neural stimulation
fhandle = @(l_CE_sol)init_muscle_force_equilib_tendon(l_CE_sol, y0(61),...
    q_sol, muscle_sol);
%l_CE_sol = fzero(fhandle, [0 y0(61)])

%l_CE_sol = fzero(fhandle, [0.04 0.05])


l_CE_sol = 0.046774811894118;
%l_CE_sol = 0.291809869606907 - 0.245;
%l_CE_sol_initial = 0.291809869606907 - 0.245;
clear fhandle

y0(62) = l_CE_sol;
% y63 is dot_l_mtc, y64 is dot_l_CE_sol, 0 at start

% TA
q_ta = 0.001;
y0(65) = l_mtc_3points(y0(15:21), y0(8:14), u_ta3, u_ta2, u_ta1);

fhandle = @(l_CE_ta)init_muscle_force_equilib(l_CE_ta, y0(65),...
    q_ta, muscle_ta);
l_CE_ta = fzero(fhandle, [0 y0(65)]);
clear fhandle

y0(66) = l_CE_ta;

% y0 69 and 70 sol length error and error dot, 0 at start
% y0 71 and 72 ta length errors

%%%

% Simulation settings
dt = 0.0002;
t_end = 0.002;
tspan = 0:dt:t_end;
options = odeset('Stats','on','RelTol',1e-6);
parameters(35) = dt;
parameters(36) = t_end;

% Simulation 
tic
[t,y] = ode45(@eomLeg, tspan, y0, options, parameters);
%[t,y] = ode15s(@eomLeg, tspan, y0, options, parameters);
toc

% gradients from integration,output values that do not need intergration
q_ta_vec =  gradient(y(:,73))/dt;
q_sol_vec = gradient(y(:,74))/dt;
F_MTC_ta =  gradient(y(:,75))/dt;
F_MTC_sol = gradient(y(:,76))/dt;
Displ_vec = gradient(y(:,77))/dt;
Tendon_strain_vec = Displ_vec/0.245;

%%%

% Post-process
%% 

global_postprocessing


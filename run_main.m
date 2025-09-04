% Multibody dynamics model of human lower limb with ANCF Achilles tendon.
% Based on the code by Erik Makkonen 2021 and Leonid Obrezkov 2022
% Aleksandr Nemov

 clc; close all; 
 clear variables;

addpath('MBS_functions','Tendon_functions')

format long
% inpuit is controlled in eomLeg.m file

% Parameters
%g = 9.81;   % m/s2
g = 0;
% Parameters are from OpenSim gait 2354 model

% Body 1 - femur 
m1 = 9.3014e4;    % kg
L11 = 0.17;     % m, distance to hip joint
L12 = 0.2258;   % m, distance to knee joint 
% Inertia tensors
Ix1_ = 0.1412;
Iy1_ = 0.0351;
Iz1_ = 0.1339;

% Body 2 - tibia and fibula
m2 = 3.7075e4;    % kg
L21 = 0.1867;   % m, distance to knee joint
L22 = 0.2433;   % m, distance to ankle joint
% Inertia tensors
Ix2_ = 0.0511;
Iy2_ = 0.0051;
Iz2_ = 0.0504;

% Body 3 - foot
% gait 2354 has 3 bodies, talus, calcaneum and toes
% that are combined into one body in this model
m3 = 1.5666;    % kg
L31 = 0.0636;   % m, length distance to ankle joint
H3 = 0.0148;    % m, height distance to ankle joint
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
t_end = 0.2;
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

%figure(2); % Animation
% Af1 = matlabFunction(body(1).A);
% Af2 = matlabFunction(body(2).A);
% Af3 = matlabFunction(body(3).A);
%for ii = 1:floor((1/dt)/200):length(t)

myVideo = VideoWriter('myVideoFile','Motion JPEG AVI'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)

figure(3); % Animation
for ii = 1:length(t)
    clf 
    mass_cm1 = [y(ii,1); y(ii,2); y(ii,3)]; 
    mass_cm2 = [y(ii,8); y(ii,9); y(ii,10)];
    mass_cm3 = [y(ii,15); y(ii,16); y(ii,17)]; 
    A1 = Af1(y(ii,4), y(ii,5), y(ii,6),y(ii,7));
    A2 = Af2(y(ii,11), y(ii,12), y(ii,13),y(ii,14));
    A3 = Af3(y(ii,18), y(ii,19), y(ii,20),y(ii,21));
      
    p1 = [0; 0; L11];
    p2 = [0; 0; -L12];
    p3 = [0; 0; L21];
    p4 = [0;  0; -L22];
    %p5 = [0; -L31; 0];
    %p6 = [0; L31; 0];
    p5 = [0; -0.1126; -0.0027];
    p6 = [0; 0.1; -0.0027];
    p1v = mass_cm1 + A1*p1;
    p2v = mass_cm1 + A1*p2;
    p3v = mass_cm2 + A2*p3;
    p4v = mass_cm2 + A2*p4;
    p5v = mass_cm3 + A3*p5;
    p6v = mass_cm3 + A3*p6;
    
    pta3 = mass_cm3 + A3*u_ta3;
    pta2 = mass_cm2 + A2*u_ta2;
    
    psol3 = mass_cm3 + A3*u_sol1;
    psol2 = mass_cm2 + A2*u_sol2;
      
   %plot3([p1v(1) p2v(1)],[p1v(2) p2v(2)],[p1v(3) p2v(3)],'r-','LineWidth',2)
   view(90,0);
   grid on
   %axis([-2 2 -1 1 -1 0])
   axis([-2 2 -0.3 0.3 -1 -0.4])
   hold on
   plot3([p3v(1) p4v(1)],[p3v(2) p4v(2)],[p3v(3) p4v(3)],'Color',[0.4 0.6 0.7],'LineWidth',5)
   plot3([p5v(1) p6v(1)],[p5v(2) p6v(2)],[p5v(3) p6v(3)],'Color',[0.5 0.5 0.5],'LineWidth',5)
   
   plot3(pta3(1),pta3(2),pta3(3), 'b*')
   plot3(pta2(1),pta2(2),pta2(3), 'b*')
   plot3([pta3(1) pta2(1)],[pta3(2) pta2(2)],[pta3(3) pta2(3)],'b-','LineWidth',1)
   
   plot3(psol3(1),psol3(2),psol3(3), 'k*')
   plot3(psol2(1),psol2(2),psol2(3), 'k*')   
   plot3([psol3(1) psol2(1)],[psol3(2) psol2(2)],[psol3(3) psol2(3)],'r-','LineWidth',1)
   
   hold off
   pause(0.05)
   frame = getframe(gcf); %get frame
   writeVideo(myVideo, frame);
end    
close(myVideo)
%% 

%plot(t,F_MTC_sol,'red',t,F_MTC_ta,'blue','LineWidth',2);
plot(t,Tendon_strain_vec*100,'red','LineWidth',2);
ax=gca;
ax.FontSize = 18;
xlim([0 0.2]);

%xlim([0 600])
%xticks(0:100:600)
%ylim([0 1.5])
%yticks(0:0.25:1.5)
xlabel('time, s','FontSize',18)
ylabel('Tendon strain, %','FontSize',18)
grid on

%legend({'Soleus','Tibialis anterior',},'FontSize',18)

figure
plot(t,q_sol_vec,'red','LineWidth',2);
% Test if total energy is constant before adding muscle forces

% Ek = zeros(1,length(t));
% Ep = zeros(1,length(t));
% E = zeros(1,length(t));
% 
% for kk = 1:length(t)
%     R31 = y(kk,3);
%     e01 = y(kk,4);
%     e11 = y(kk,5);
%     e21 = y(kk,6);
%     e31 = y(kk,7);
%     R32 = y(kk,10);
%     e02 = y(kk,11);
%     e12 = y(kk,12);
%     e22 = y(kk,13);
%     e32 = y(kk,14);
%     R33 = y(kk,17);
%     e03 = y(kk,18);
%     e13 = y(kk,19);
%     e23 = y(kk,20);
%     e33 = y(kk,21);
%     
%     Ek(kk) = 0.5*y(kk,22:42)*Mf3(Ix1_,Ix2_,Ix3_,Iy1_,Iy2_,Iy3_,...
%         Iz1_,Iz2_,Iz3_,e01,e02,e03,e11,e12,e13,e21,e22,e23,e31,...
%         e32,e33,m1,m2,m3)*y(kk,22:42).';            % 2.85
%     Ep(kk) = m1*g*(R31) + m2*g*(R32) + m3*g*R33;        % Ep = mgh
%     E(kk) = Ek(kk) + Ep(kk);
% end
% figure(1);
% plot(t,Ek, t,Ep, t,E)
% xlabel('Time [s]')
% ylabel('System energy [J]')
% legend('Kinetic', 'Potential', 'Total')
% grid on

function [dydt] = eomLeg(t,y,parameters)
% solve equations of motions

R11 = y(1);
R21 = y(2);
R31 = y(3);
e01 = y(4);
e11 = y(5);
e21 = y(6);
e31 = y(7);

R12 = y(8);
R22 = y(9);
R32 = y(10);
e02 = y(11);
e12 = y(12);
e22 = y(13);
e32 = y(14);

R13 = y(15);
R23 = y(16);
R33 = y(17);
e03 = y(18);
e13 = y(19);
e23 = y(20);
e33 = y(21);

dR11 = y(22);
dR21 = y(23);
dR31 = y(24);
de01 = y(25);
de11 = y(26);
de21 = y(27);
de31 = y(28);

dR12 = y(29);
dR22 = y(30);
dR32 = y(31);
de02 = y(32);
de12 = y(33);
de22 = y(34);
de32 = y(35);

dR13 = y(36);
dR23 = y(37);
dR33 = y(38);
de03 = y(39);
de13 = y(40);
de23 = y(41);
de33 = y(42);

g = parameters(1);

m1 = parameters(2);
L11 = parameters(3);
L12 = parameters(4);
Ix1_ = parameters(5);
Iy1_ = parameters(6);
Iz1_ = parameters(7);

m2 = parameters(8);
L21 = parameters(9);
L22 = parameters(10);
Ix2_ = parameters(11);
Iy2_ = parameters(12);
Iz2_ = parameters(13);

m3 = parameters(14);
L31 = parameters(15);
H_3 = parameters(16);
Ix3_ = parameters(17);
Iy3_ = parameters(18);
Iz3_ = parameters(19);

dt = parameters(35);
t_end = parameters(36);

body2_cm = [y(8); y(9); y(10)];
body3_cm = [y(15); y(16); y(17)];

A2 = Af2(y(11), y(12), y(13), y(14));
A3 = Af3(y(18), y(19), y(20), y(21));

u_p1 = [0 0.10 0].';
sim_p1 = body3_cm + A3*u_p1;

q_min = 0.001;

% index for current time

index = round(t*1/dt+1);

% length of the vector must be the same as time
tspan = 0:dt:t_end;


% Display current time
fprintf('Current time: %.4f s (index %d of %d)\n', t, index, length(tspan));


input = 5;
% 1 is stable 
% 2 is moving up 
% 3 is moving down
% 4 is moving up and down, MORE WORK NEEDED
% 5 direct input of a q pulse

if input == 1
    input_points = [-0.8406 -0.8406];
    input_ptime = linspace(0, t_end, length(input_points));
    input_curve = pchip(input_ptime, input_points, tspan);
elseif input == 2
    %input_points = [-0.8406 -0.8406 -0.8306 -0.8306 -0.8206 -0.8206];
    input_points = [-0.8406 -0.8206 -0.8206];
    input_ptime = linspace(0, t_end, length(input_points));
    input_curve = pchip(input_ptime, input_points, tspan);
elseif input == 3
    input_points = [-0.8406 -0.8706 -0.8706];
    %input_points = [-0.8406 -0.8706 -0.8706 -0.8406 -0.8406];
    input_ptime = linspace(0, t_end, length(input_points));
    input_curve = pchip(input_ptime, input_points, tspan);
    %input_points = [-0.8406 -0.8406 -0.8406 -0.8706 -0.8706 -0.8706];
    %input_ptime = linspace(0, t_end, length(input_points));
    %input_curve = pchip(input_ptime, input_points, tspan);
elseif input == 4
    input_points = [-0.8406 -0.8406 -0.8206 -0.8206 -0.8606 -0.8606 ...
        -0.8406 -0.8406 -0.8606 -0.8606 -0.8306 -0.8306 ...
        -0.8406 -0.8406];
    input_ptime = linspace(0, t_end, length(input_points));
    input_curve = pchip(input_ptime, input_points, tspan);
elseif input == 5
    %input_points = [0.001 0.001 1 1 0.001 0.001];
    
    %input_points = [0.001 0.2 0.2 0.001 0.001];
    input_points = [0.001 0.2];

    %input_points = [0.001 1 1 0.001 0.001];
    %input_points = [0.001 0.001];
    input_ptime = linspace(0,t_end,length(input_points));
    input_curve = pchip(input_ptime,input_points,tspan);   
end

% muscle components
u_sol1 = parameters(20:22).';
u_sol2 = parameters(23:25).';
muscle_sol = muscle_param_sol();

l_MTC_sol_int = y(61);
l_MTC_sol = l_mtc_2points(y(15:21), y(8:14), u_sol1, u_sol2);
Diff = l_MTC_sol - l_MTC_sol_int;
l_CE_sol = y(62);
%dot_l_MTC_sol = 0;
dot_l_MTC_sol = dot_l_mtc_2points(y(15:21), y(8:14), y(36:42), y(29:35), u_sol1, u_sol2);

u_ta1 = parameters(26:28).';
u_ta2 = parameters(29:31).';
u_ta3 = parameters(32:34).';
muscle_ta = muscle_param_ta();

l_MTC_ta = l_mtc_3points(y(15:21), y(8:14), u_ta3, u_ta2, u_ta1);
l_CE_ta = y(66);
%dot_l_MTC_ta = 0;
dot_l_MTC_ta = dot_l_mtc_3points(y(15:21), y(8:14), y(36:42), y(29:35), u_ta3, u_ta2, u_ta1);

% Common PID CONTROLLER
%period_sol = 10;
%kp_sol = 40;                    % proportional gain
%ki_sol = period_sol*0.5;        % integrate gain 50% of period 
%kd_sol = period_sol*0.125*2;    % derivative gain 12.5%, 25%

% period_ta = 40;
% kp_ta = 30;                     % proportional gain
% ki_ta = period_ta*0.5;          % integrate gain 50% of period 
% kd_ta = period_ta*0.125*2;      % derivative gain 12.5%, 25%

% %sim_p1(3)
% %input_curve(index)
% error_p = sim_p1(3) - input_curve(index);
% error_i = y(69);
% error_d = (error_p - y(70))/dt;

% activity

q_sol = input_curve(index);
q_ta = q_min;

% if error_p > 0
%     % reduce proportional gain near target
%     if error_p < 5e-4
%         kp_sol = kp_sol*abs(error_d)/5e-4;
%     end
%     q_sol = min(1,kp_sol*error_p + kd_sol*error_d + ki_sol*error_i)
%     %if (q_sol < q_min)  || (error_p < 0)
%     if (q_sol < q_min)  
%         q_sol = q_min;
%     end
%     q_ta = q_min;
% else 
%     % reduce proportional gain near target
%     if abs(error_p) < 5e-4
%         kp_ta = kp_ta*abs(error_d)/5e-4;
%     end
%     q_ta = min(1,kp_ta*(-error_p) + kd_ta*(-error_d) + ki_ta*(-error_i));
%     %if (q_ta < q_min) || (error_p < -1e-4)
%     if (q_ta < q_min) 
%         q_ta = q_min;
%     end
%     q_sol = q_min;
% end

% PID CONTROLLER (SOL)
%period_sol = 10;
%kp_sol = 40;                    % proportional gain
%ki_sol = period_sol*0.5;        % integrate gain 50% of period 
%kd_sol = period_sol*0.125*2;    % derivative gain 12.5%, 25%

%error_p_sol = sim_p1(3) - input_curve(index);
%error_i_sol = y(69);
%error_d_sol = (error_p_sol - y(70))/dt;

% reduce proportional gain near target
%if error_p_sol < 5e-4
%    kp_sol = kp_sol*error_d_sol/5e-4;
%end

% activity
%q_sol = min(1,kp_sol*error_p_sol + kd_sol*error_d_sol + ki_sol*error_i_sol);
%if (q_sol < q_min)  || (error_p_sol < 0)
%    q_sol = q_min;
%end
%q_sol = q_min;
%q_sol = 0.2*t/t_end

% PID CONTROLLER (TA)
%period_ta = 40;
%kp_ta = 30;                     % proportional gain
%ki_ta = period_ta*0.5;          % integrate gain 50% of period 
%kd_ta = period_ta*0.125*2;      % derivative gain 12.5%, 25%

%error_p_ta = -error_p_sol;
%error_i_ta = -error_i_sol;
%error_d_ta = -error_d_sol;

%if error_p_ta < 5e-4
%    kp_ta = kp_ta*error_d_ta/5e-4;
%end

%q_ta = min(1,kp_ta*error_p_ta + kd_ta*error_d_ta + ki_ta*error_i_ta);
%if (q_ta < q_min) || (error_p_ta < -1e-4)
%    q_ta = q_min;
%end
%q_ta = 0.1*t/t_end
%q_ta = q_min


% foot velocity y 37 and 38
% if y(38) < -0.025
%     y(38)
%     q_ta = min(1,q_ta + 1.5*(-y(38)-0.025));
% end
% 
% if y(38) > 0.025
%     y(38)
%     q_sol = min(1,q_sol + 1.5*(y(38)-0.025));
% end

% MUSCLE FORCES
% SOL force


%[F_MTC_sol, dot_l_CE_sol, F_elements_sol] = mtc_model_matlab(l_CE_sol,...
[F_MTC_sol, dot_l_CE_sol, F_elements_sol, Displ] = mtc_model_w_tendon(l_CE_sol,...
    l_MTC_sol, dot_l_MTC_sol, q_sol, muscle_sol);

%F_MTC_sol = 0.5305;
%dot_l_CE_sol = 1.1166;
%F_elements_sol = [27.1644
%         0
%    0.9108
%    0.5305
%  -26.6339];
% force direction
p1_sol = body3_cm + A3*u_sol1;
p2_sol = body2_cm + A2*u_sol2;
dp_sol = p2_sol - p1_sol;

Fx_sol = 0;
Fy_sol = F_MTC_sol*dp_sol(2)/sqrt(sum(dp_sol.^2));
Fz_sol = F_MTC_sol*dp_sol(3)/sqrt(sum(dp_sol.^2));

% TA force

[F_MTC_ta, dot_l_CE_ta, F_elements_ta] = mtc_model_matlab(l_CE_ta, ...
    l_MTC_ta, dot_l_MTC_ta, q_ta, muscle_ta);
%F_MTC_sol
%F_MTC_ta

% force direction
% body 2
p1_ta = body2_cm + A2*u_ta1;
p2_ta = body2_cm + A2*u_ta2;
dp1_ta = p1_ta - p2_ta;
Fx_ta2 = 0;
Fy_ta2 = F_MTC_ta*dp1_ta(2)/sqrt(sum(dp1_ta.^2));
Fz_ta2 = F_MTC_ta*dp1_ta(3)/sqrt(sum(dp1_ta.^2));
% body 3
p3_ta = body3_cm + A3*u_ta3;
dp2_ta = p2_ta - p3_ta;
Fx_ta3 = 0;
Fy_ta3 = F_MTC_ta*dp2_ta(2)/sqrt(sum(dp2_ta.^2));
Fz_ta3 = F_MTC_ta*dp2_ta(3)/sqrt(sum(dp2_ta.^2));

% System matrices
sysM = SysP3Mf3(H_3,Ix1_,Ix2_,Ix3_,Iy1_,Iy2_,Iy3_,Iz1_,Iz2_,Iz3_,...
    L11,L12,L21,L22,L31,R11,R12,R13,R21,R22,R23,R31,R32,R33,...
    e01,e02,e03,e11,e12,e13,e21,e22,e23,e31,e32,e33,m1,m2,m3);

sysF = SysP3Ff3(Fx_ta2,Fx_ta3,Fx_sol,Fy_ta2,Fy_ta3,Fy_sol,Fz_ta2,Fz_ta3,...
    Fz_sol,H_3,Ix1_,Ix2_,Ix3_,Iy1_,Iy2_,Iy3_,...
    Iz1_,Iz2_,Iz3_,L11,L12,L21,L22,L31,R11,R12,R13,R21,R22,R23,R31,...
    R32,R33,dR11,dR12,dR13,dR21,dR22,dR23,dR31,dR32,dR33,de01,de02,...
    de03,de11,de12,de13,de21,de22,de23,de31,de32,de33,e01,e02,e03,...
    e11,e12,e13,e21,e22,e23,e31,e32,e33,g,m1,m2,m3);

% EOM
EOM = sysM\sysF;

dydt(1:21) = y(22:42);
dydt(22:60) = EOM;

dydt(61:62) = y(63:64);
dydt(63) = dot_l_MTC_sol;
dydt(64) = dot_l_CE_sol;
dydt(65:66) = y(67:68);
dydt(67) = 0;
dydt(68) = dot_l_CE_ta;
dydt(69) = y(70);
%dydt(70) = error_d;
dydt(70) = 0;
dydt(71) = y(72);
%dydt(72) = error_d;
dydt(72) = 0;
dydt(73) = q_ta;
dydt(74) = q_sol;
dydt(75) = F_MTC_ta;
dydt(76) = F_MTC_sol;
dydt(77) = Displ;

dydt = dydt(:);

end

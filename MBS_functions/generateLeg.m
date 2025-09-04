function body = generateLeg()
% Generating system matrices using Euler parameters

% symbolic variables
syms R11 R21 R31 e01 e11 e21 e31
syms dR11 dR21 dR31 de01 de11 de21 de31
syms R12 R22 R32 e02 e12 e22 e32            
syms dR12 dR22 dR32 de02 de12 de22 de32
syms R13 R23 R33 e03 e13 e23 e33
syms dR13 dR23 dR33 de03 de13 de23 de33
syms g m1 m2 m3 L11 L12 L21 L22 L31 H_3
syms Fx_sol Fy_sol Fz_sol
syms Fx_ta2 Fy_ta2 Fz_ta2
syms Fx_ta3 Fy_ta3 Fz_ta3

syms Ix1_ Iy1_ Iz1_     
syms Ix2_ Iy2_ Iz2_
syms Ix3_ Iy3_ Iz3_

nBodies = 3;

% properites 
body(1).m = m1;
body(1).L = L11+L12;

body(2).m = m2;
body(2).L = L21+L22;

body(3).m = m3;
body(3).L = 2*L31;
body(3).H = H_3;

body(1).Izz_ = Iz1_;
body(1).Ixx_ = Ix1_;
body(1).Iyy_ = Iy1_;

body(2).Izz_ = Iz2_;
body(2).Ixx_ = Ix2_;
body(2).Iyy_ = Iy2_;

body(3).Izz_ = Iz3_;
body(3).Ixx_ = Ix3_;
body(3).Iyy_ = Iy3_;


% coordinate vectors
body(1).q = [R11, R21, R31, e01 e11 e21 e31].';
body(2).q = [R12, R22, R32, e02 e12 e22 e32].';
body(3).q = [R13, R23, R33, e03 e13 e23 e33].';

body(1).dq = [dR11, dR21, dR31, de01 de11 de21 de31].';
body(2).dq = [dR12, dR22, dR32, de02 de12 de22 de32].';
body(3).dq = [dR13, dR23, dR33, de03 de13 de23 de33].';

body(1).e = [e11; e21; e31];
body(1).de = [de11; de21; de31];
body(1).e0 = e01;
body(1).de0 = de01;

body(2).e = [e12; e22; e32];
body(2).de = [de12; de22; de32];
body(2).e0 = e02;
body(2).de0 = de02;

body(3).e = [e13; e23; e33];
body(3).de = [de13; de23; de33];
body(3).e0 = e03;
body(3).de0 = de03;

q = [body(1).q; body(2).q; body(3).q];
dq = [body(1).dq; body(2).dq; body(3).dq];

% Rotation matrices
for ii = 1:nBodies
    body(ii).p = [body(ii).e0; body(ii).e];
    body(ii).dp = [body(ii).de0; body(ii).de];
    
    body(ii).A = (2*body(ii).e0^2 - 1)*eye(3) + ...
        2*(body(ii).e*body(ii).e.' + body(ii).e0*skew(body(ii).e));

    body(ii).E = [-body(ii).e, skew(body(ii).e) + body(ii).e0*eye(3)];  
    body(ii).E_ = [-body(ii).e, -skew(body(ii).e) + body(ii).e0*eye(3)];
    
    body(ii).G = 2*body(ii).E;
    body(ii).G_ = 2*body(ii).E_;
end

% Angular velocity of the system 
for ii = 1:nBodies
    body(ii).w_ = body(ii).G_*body(ii).dp;
end

% Mass matrix
for ii = 1:nBodies
    body(ii).Itt_ = [
        body(ii).Ixx_, 0, 0;
        0, body(ii).Iyy_, 0;
        0, 0, body(ii).Izz_];

    body(ii).Mrr = body(ii).m * eye(3);

    body(ii).Mtt = body(ii).G_.' * body(ii).Itt_ * body(ii).G_;

    body(ii).M(1:3,1:3) = body(ii).Mrr;
    body(ii).M(4:7,4:7) = body(ii).Mtt;
end

% Constraints
% Euler parameter constraints
for ii = 1:nBodies
    body(ii).c = body(ii).dp.'*body(ii).dp;
    body(ii).P = [zeros(1,3), body(ii).p.'];
end

% Revolute joint y1_ = -L/2 , X rotation allowed
C(1:3) = body(1).q(1:3) + body(1).A*[0; 0; L11];

vec21 = body(1).q(1:3) + body(1).A*[0; 100; 0];  % y direction 1 body
vec31 = body(1).q(1:3) + body(1).A*[0; 0; 100];  % z direction 1 body

cp1 = skew(vec21)*vec31;  
  
C(4) = simplify(cp1(3));
C(5) = simplify(cp1(2));

% Revolute joint connecting bodies 1 and 2
C(6:8) = body(1).q(1:3) + body(1).A*[0; 0; -L12] - ...
    (body(2).q(1:3) + body(2).A*[0; 0; L21]);

vec22 = body(2).q(1:3) + body(2).A*[0; 100; 0];  % y direction 2 body
vec32 = body(2).q(1:3) + body(2).A*[0; 0; 100];  % z direction 2 body

cp2 = skew(vec22)*vec32;  
C(9) = simplify(cp2(3));
C(10) = simplify(cp2(2));

% Revolute joint connecting bodies 2 and 3
C(11:13) = body(2).q(1:3) + body(2).A*[0; 0; -L22] - ...
    (body(3).q(1:3) + body(3).A*[0; -L31; H_3]);

vec23 = body(3).q(1:3) + body(3).A*[0; 100; 0];  % y direction 3 body
vec33 = body(3).q(1:3) + body(3).A*[0; 0; 100];  % z direction 3 body

cp3 = skew(vec23)*vec33;  
C(14) = simplify(cp3(3));
C(15) = simplify(cp3(2));

C = C(:);
Cq = jacobian(C, q);
Qc = -jacobian(Cq*dq,q)*dq;

% Appy Baumgarten stabilisation for constraints
alpha = 1/1e-3;      % assume dt = 1e-3
beta = sqrt(2)/1e-3;   
Qc = Qc - (2*alpha*Cq*dq + beta^2*C);


% Vector of external forces (muscles)

% m. soleus force 
body(2).F_sol = [-Fx_sol; -Fy_sol; -Fz_sol]; 
body(3).F_sol = [Fx_sol; Fy_sol; Fz_sol];
% force location 
u_sol21 = [0 -0.002 0.0337];  % tibia
u_sol31 = [0 -0.1126 0.0038]; % heel

% m. TA force
body(2).F_ta = [Fx_ta2; Fy_ta2; Fz_ta2];
body(3).F_ta = [Fx_ta3; Fy_ta3; Fz_ta3];
% force location
u_ta21 = [0 0.018 0.0247];      % tibia
u_ta31 = [0 0.0044 -0.0092];    % foot

body(1).Q_trans_ = [0; 0; -body(1).m*g];
body(2).Q_trans_ = [0; 0; -body(2).m*g] + body(2).F_sol + body(2).F_ta;
body(3).Q_trans_ = [0; 0; -body(3).m*g] + body(3).F_sol + body(3).F_ta;

body(1).Qe = [body(1).Q_trans_; zeros(4,1)];

body(2).Qe = [body(2).Q_trans_; ...
    body(2).G_.'*skew(u_sol21)*body(2).A.'*body(2).F_sol + ...
    body(2).G_.'*skew(u_ta21)*body(2).A.'*body(2).F_ta];

body(3).Qe = [body(3).Q_trans_; ...
    body(3).G_.'*skew(u_sol31)*body(3).A.'*body(3).F_sol + ...
    body(3).G_.'*skew(u_ta31)*body(3).A.'*body(3).F_ta];    

for ii = 1:nBodies
    % Quadratic velocity vector
    body(ii).dE_ = diff(body(ii).E_, body(ii).q(4))*body(ii).dq(4)...
        + diff(body(ii).E_, body(ii).q(5))*body(ii).dq(5)...
        + diff(body(ii).E_, body(ii).q(6))*body(ii).dq(6)...
        + diff(body(ii).E_, body(ii).q(7))*body(ii).dq(7); % E_dash_dot

    body(ii).Qvt = 2*2*body(ii).dE_.'*body(ii).Itt_*body(ii).w_;

    body(ii).Qv(1:3,1) = sym(zeros(3,1));
    body(ii).Qv(4:7,1) = body(ii).Qvt;
end

% System matrices
M = sym(zeros(7*nBodies));
P = sym(zeros(nBodies, 7*nBodies));

for ii = 1:nBodies
    M((ii-1)*7+1:ii*7,(ii-1)*7+1:ii*7) = body(ii).M; % mass matrix
    P(ii,(ii-1)*7+1:ii*7) = body(ii).P;   % Euler paremeter constraints
    
    Qe((ii-1)*7+1:ii*7) = body(ii).Qe; % external forces
    Qv((ii-1)*7+1:ii*7) = body(ii).Qv; % quadratic forces    
    c(ii) = body(ii).c;                % Constraint forces related to P
end
     
sysM = simplify([M, P.', Cq.';
    P, zeros(nBodies), zeros(nBodies,length(C));
    Cq, zeros(length(C),nBodies), zeros(length(C))]);

sysF = simplify([Qe.' - Qv.';
    -c.';
    Qc]);

% Generate function files
matlabFunction(sysM, 'file', 'MBS_functions/SysP3Mf3');
matlabFunction(sysF, 'file', 'MBS_functions/SysP3Ff3');
matlabFunction(M,'file','MBS_functions/Mf3');
matlabFunction(C,'file','MBS_functions/Cf3'); 
matlabFunction(Cq, 'file', 'MBS_functions/Cqf3');

matlabFunction(body(1).A, 'file', 'MBS_functions/Af1');
matlabFunction(body(2).A, 'file', 'MBS_functions/Af2');
matlabFunction(body(3).A, 'file', 'MBS_functions/Af3');

end
% Multibody dynamics model of human lower limb with ANCF Achilles tendon.
% Based on the code by Erik Makkonen 2021 and Leonid Obrezkov 2022
% Aleksandr Nemov

function run_main

 clc; close all; 
 clear variables;

addpath('MBS_functions','Tendon_functions')
addpath("Postprocessing");

format long
% inpuit is controlled in eomLeg.m file

% Parameters
%g = 9.81;   % m/s2
g = 0;
% Parameters are from LaiArnold2017 model

% Body 1 - femur 
m1 = 9.3014e4;    % kg
L11 = 0.17;     % m, distance to hip joint
L12 = 0.2396;   % m, distance to knee joint 

% Inertia tensors
Ix1_ = 0.1412;
Iy1_ = 0.0351;
Iz1_ = 0.1339;

% Body 2 - tibia and fibula
% Parameters from 
m2 = 3.7075e4;    % kg
L21 = 0.1867;   % m, distance to knee joint
L22 = 0.2133;   % m, distance to ankle joint
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
L31 = 0.04;
%H3 = 0.0148;    % m, height distance to ankle joint
H3 = 0.03;    % m, height distance to ankle joint
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

l_CE_sol = 0.037165625152305;
%l_CE_sol = 0.046774811894118;
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
dt = 0.001;
t_end = 0.6;
tspan = 0:dt:t_end;
options = odeset('Stats','on','RelTol',1e-6);
parameters(35) = dt;
parameters(36) = t_end;

% Add storage for tendon deformation data
Body1_q_values = [];
Body2_q_values = [];
Body3_q_values = [];


% Simulation 
tic
[t,y] = ode45(@eomLegWrapper, tspan, y0, options);
%[t,y] = ode45(@eomLeg, tspan, y0, options, parameters);
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

% Get full Body1, Body2, Body3 structures for visualization
[~, Body1, Body2, Body3] = eomLeg(t(end), y(end,:)', parameters);

% Save entire workspace
save('AchillesWorkspace.mat');


% Post-process
%% 

% Get screen size
screenSize = get(0, 'ScreenSize');

% Create a new figure and set its position to full screen
figure('Position', [0 0 screenSize(3) screenSize(4)]);

% Import the STL files
verticesShin_initial= stlread("tibia_r_from_OpenSim_smooth.stl");
verticesFemur_initial = stlread("femur_r_from_OpenSim.stl");
verticesFoot_initial = stlread("foot_from_OpenSim_fix_smooth.stl");


% Define custom scale factors
scaleShin = [1, 1.1, 1];  % scale X, Y, Z

% Apply scaling to tibia vertices
scaledVerticesShin = verticesShin_initial.Points .* scaleShin;

% Apply rotation around X-axis by -90 degrees
%R = eye(3);
%R = [1 0 0; 0 0 1; 0 -1 0];
R = [0 0 1; 1 0 0; 0 1 0];
verticesShinPoints = (R * 1000* scaledVerticesShin')';
verticesFemurPoints = (R * 1000 * verticesFemur_initial.Points')';
verticesFootPoints = (R * 1000* verticesFoot_initial.Points')';

% Extract points and connectivity lists

% Center mass coordinates in local CS
femur_local_CM = 1000*[ 1.5943e-05; 4.97468e-06; -0.000206997];
shin_local_CM = 1000*[6.9501e-06; 5.79781e-06; -0.000196268*1.2];
foot_local_CM = 1000*[4.01932e-06; 8.01158e-05; 1.59065e-05];

shinCOM=shin_local_CM'*1000; % point from the ankle joint to COM
pointsShin = verticesShinPoints-shinCOM;
connectivityShin = verticesShin_initial.ConnectivityList;

%femurCOM=[-0.633289;-0.802438;-1.08622].';
femurCOM=femur_local_CM'*1000;
pointsFemur = verticesFemurPoints-femurCOM;
connectivityFemur = verticesFemur_initial.ConnectivityList;

%footCOM=[52.7643;-45.9727;2.72611].';
footCOM=foot_local_CM'*1000;
pointsFoot = verticesFootPoints-footCOM;
connectivityFoot = verticesFoot_initial.ConnectivityList;


i=1;


%% Prepare video writers for writing animations
writerFull = VideoWriter('Plantar_flextion_full.mp4', 'MPEG-4');
writerZoom = VideoWriter('Plantar_flextion_zoom.mp4', 'MPEG-4');
writerFull.FrameRate = 10;
writerZoom.FrameRate = 10;
open(writerFull);
open(writerZoom);



for ii = 1:size(y,1)
%for ii = 1:1
    if isnan(y(ii,1))
        break;
    end
    title('Transformed Geometry');
    axis equal;
    grid on;  


    % Calculate dynamic transformations
    translationFemur = y(ii,1:3)*1000;
    translationShin = y(ii,1+7:3+7)*1000;
    translationFoot = y(ii,1+14:3+14)*1000;

    RFemur=Rotationmatrix_euler_parameter(y(ii,4),y(ii,5),y(ii,6),y(ii,7));
    RShin=Rotationmatrix_euler_parameter(y(ii,4+7),y(ii,5+7),y(ii,6+7),y(ii,7+7));
    RFoot=Rotationmatrix_euler_parameter(y(ii,4+2*7),y(ii,5+2*7),y(ii,6+2*7),y(ii,7+2*7));

    % Transform shin vertices
    transformedVerticesFemur = (RFemur * pointsFemur')';
    transformedVerticesFemur = transformedVerticesFemur + repmat(translationFemur, length(transformedVerticesFemur), 1);

    % Transform Foot vertices
    transformedVerticesShin = (RShin * pointsShin')';
    transformedVerticesShin = transformedVerticesShin + repmat(translationShin, length(transformedVerticesShin), 1);

   
    % Transform Foot vertices
    transformedVerticesFoot = (RFoot * pointsFoot')';
    transformedVerticesFoot = transformedVerticesFoot + repmat(translationFoot, length(transformedVerticesFoot), 1);

    % Create updated geometry for foot
    modelFemur = createpde();
    geometryFromMesh(modelFemur, transformedVerticesFemur', connectivityFemur');
    pdegplot(modelFemur);
    set(findobj(gca,'Type','Quiver'),'Visible','off');
    set(findall(gca,'Layer','Middle'),'Visible','Off');
    %view(90,0);

    hold on;

    % Create updated geometry for shin
    modelShin = createpde();
    geometryFromMesh(modelShin, transformedVerticesShin', connectivityShin');
    pdegplot(modelShin);
    set(findobj(gca,'Type','Quiver'),'Visible','off');
    set(findall(gca,'Layer','Middle'),'Visible','Off');

    view(90,0);
    %axis([-100 100 -300 200 -700 100])  

    hold on;
    % Create updated geometry for foot
    modelFoot = createpde();
    geometryFromMesh(modelFoot, transformedVerticesFoot', connectivityFoot');
    pdegplot(modelFoot);
    set(findobj(gca,'Type','Quiver'),'Visible','off');
    set(findall(gca,'Layer','Middle'),'Visible','Off');
    view(90,0);
    %axis([-100 100 -300 200 -700 100])  
    %MATLAB Script to Convert Figure to Image

    xlabel('X');
    ylabel('Y');
    zlabel('Z');


    % Adding visualization of the deformed tendon

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
     p5 = [0; -L31; 0];
     p6 = [0; L31; 0];
     % p5 = [0; -0.1126; -0.0027];
     % p6 = [0; 0.1; -0.0027];
     p1v = mass_cm1 + A1*p1;
     p2v = mass_cm1 + A1*p2;
     p3v = mass_cm2 + A2*p3;
     p4v = mass_cm2 + A2*p4;
     p5v = mass_cm3 + A3*p5;
     p6v = mass_cm3 + A3*p6;
%     pta3 = mass_cm3 + A3*u_ta3;
%     pta2 = mass_cm2 + A2*u_ta2;
    u_sol1 = [0; -0.086; 0.003];
     psol3 = mass_cm3 + A3*u_sol1;
     psol2 = mass_cm2 + A2*u_sol2;
     tendon_visualization(Body1,Body1_q_values(ii,:)'*1000,'cyan',true,(psol2-psol3)*1000,psol3*1000);
     tendon_visualization(Body2,Body2_q_values(ii,:)'*1000,'red',true,(psol2-psol3)*1000,psol3*1000);
     tendon_visualization(Body3,Body3_q_values(ii,:)'*1000,'blue',true,(psol2-psol3)*1000,psol3*1000);
     plot3([psol3(1)*1000 psol2(1)*1000],[psol3(2)*1000 psol2(2)*1000],[psol3(3)*1000 psol2(3)*1000],'r-','LineWidth',1)

    hold off;
    %% Capture full view frame
    frameFull = getframe(gcf);
    writeVideo(writerFull, frameFull);

    %% Apply zoom for second video
    zoomCenter = (psol3 + (psol2-psol3)/5) * 1000;  % midpoint in mm
    %zoomCenter = psol3 * 1000;  % midpoint in mm
    zoomRange = 50;  % mm
    % xlim([zoomCenter(1)-zoomRange, zoomCenter(1)+zoomRange]);
    % ylim([zoomCenter(2)-zoomRange, zoomCenter(2)+zoomRange]);
    % zlim([zoomCenter(3)-zoomRange, zoomCenter(3)+zoomRange]);

    camtarget(zoomCenter);  % Set the camera target to tendon midpoint
    camzoom(10);             % Zoom in by factor of 2 (adjust as needed)


    frameZoom = getframe(gcf);
    writeVideo(writerZoom, frameZoom);

    % F(i) = getframe(gcf) ;
    i=i+1;
    drawnow
end


%% Close both writers
close(writerFull);
close(writerZoom);


% writerObj = VideoWriter('myVideo.mp4', 'MPEG-4');
% writerObj.FrameRate = 10;
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%      % convert the image to a frame
%      frame = F(i) ;    
%      writeVideo(writerObj, frame);
%  end
% % close the writer object
% close(writerObj);


%% Nested function enabling storing deformation data from all simulation steps

function dydt = eomLegWrapper(t, y)
    [dydt, Body1, Body2, Body3] = eomLeg(t, y, parameters);
    Body1_q_values = [Body1_q_values; Body1.q(:)']; 
    Body2_q_values = [Body2_q_values; Body2.q(:)']; 
    Body3_q_values = [Body3_q_values; Body3.q(:)'];  
end

end
% Get screen size
screenSize = get(0, 'ScreenSize');

% Create a new figure and set its position to full screen
figure('Position', [0 0 screenSize(3) screenSize(4)]);

% Import the STL files
verticesShin_initial= stlread("Human_tibia_and_fibula_3.STL");
verticesFemur_initial = stlread("femur_r.stl");
verticesFoot_initial = stlread("calcn_3.STL");

% Apply rotation around X-axis by -90 degrees
%R = eye(3);
%R = [1 0 0; 0 0 1; 0 -1 0];
R = [0 0 1; 1 0 0; 0 1 0];
verticesShinPoints = (R * verticesShin_initial.Points')';
verticesFemurPoints = (R * verticesFemur_initial.Points')';
verticesFootPoints = (R * verticesFoot_initial.Points')';

% Extract points and connectivity lists

% Center mass coordinates in local CS
femur_local_CM = [-0.00108622; -0.000633289; -0.000802438];
shin_local_CM = [0.000499889; 0.00312171; 0.210276];
foot_local_CM = [0.00272611; 0.0527643; -0.0459727];
%foot_local_CM = [0; 0; 0];

% Joints location in local CS

%femurCS_to_shin = [3.71941/1000; 210.576/1000; -237.488/1000];
%shinCS_to_femur = [2.05448/1000; 565.884/1000; 369.835/1000];
%shinCS_to_foot = [-1.57871/1000; 598.252/1000; -7.11398/1000];
%footCS_to_shin = [-2.70888/1000; 824.845/1000; 63.5783/1000];


%shinCOM=[3.12171;210.276;0.49989].'; % point from the ankle joint to COM
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



y = zeros(1,size(y0,1));
y = y0';

i=1;

% for ii = 1:10:size(y,2)
for ii = 1:1
    if isnan(y(1,ii))
        break;
    end
    title('Transformed Geometry');
    axis equal;
    grid on;
    %view(2);


    % Calculate dynamic transformations
    translationFemur = y(ii,1:3)*1000
    % translationFemur = [0 0 0];
    translationShin = y(ii,1+7:3+7)*1000
    % translationShin = [0 0 0];
    translationFoot = y(ii,1+14:3+14)*1000
    %translationFoot = [0 0 0];

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
    view(90,0);

    % axis([-100 100 -300 200 -700 100])  
    % MATLAB Script to Convert Figure to Image

    hold on;
    % 
    % Create updated geometry for shin
    modelShin = createpde();
    geometryFromMesh(modelShin, transformedVerticesShin', connectivityShin');
    pdegplot(modelShin);
    set(findobj(gca,'Type','Quiver'),'Visible','off');
    set(findall(gca,'Layer','Middle'),'Visible','Off');

    view(90,0);
    % % axis([-100 100 -300 200 -700 100])  
    % 
    hold on;
    % Create updated geometry for foot
    modelFoot = createpde();
    geometryFromMesh(modelFoot, transformedVerticesFoot', connectivityFoot');
    pdegplot(modelFoot);
    set(findobj(gca,'Type','Quiver'),'Visible','off');
    set(findall(gca,'Layer','Middle'),'Visible','Off');
    %view(90,0);
    % axis([-100 100 -300 200 -700 100])  
    %MATLAB Script to Convert Figure to Image

    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold off;
    F(i) = getframe(gcf) ;
    i=i+1;
    drawnow
end

writerObj = VideoWriter('myVideo.mp4', 'MPEG-4');
writerObj.FrameRate = 10;
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);



%%%%% Old postprocessing script below %%%%%%

%figure(2); % Animation
% Af1 = matlabFunction(body(1).A);
% Af2 = matlabFunction(body(2).A);
% Af3 = matlabFunction(body(3).A);
%for ii = 1:floor((1/dt)/200):length(t)

% myVideo = VideoWriter('myVideoFile','Motion JPEG AVI'); %open video file
% myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
% open(myVideo)
% 
% figure(3); % Animation
% for ii = 1:length(t)
%     clf 
%     mass_cm1 = [y(ii,1); y(ii,2); y(ii,3)]; 
%     mass_cm2 = [y(ii,8); y(ii,9); y(ii,10)];
%     mass_cm3 = [y(ii,15); y(ii,16); y(ii,17)]; 
%     A1 = Af1(y(ii,4), y(ii,5), y(ii,6),y(ii,7));
%     A2 = Af2(y(ii,11), y(ii,12), y(ii,13),y(ii,14));
%     A3 = Af3(y(ii,18), y(ii,19), y(ii,20),y(ii,21));
% 
%     p1 = [0; 0; L11];
%     p2 = [0; 0; -L12];
%     p3 = [0; 0; L21];
%     p4 = [0;  0; -L22];
%     %p5 = [0; -L31; 0];
%     %p6 = [0; L31; 0];
%     p5 = [0; -0.1126; -0.0027];
%     p6 = [0; 0.1; -0.0027];
%     p1v = mass_cm1 + A1*p1;
%     p2v = mass_cm1 + A1*p2;
%     p3v = mass_cm2 + A2*p3;
%     p4v = mass_cm2 + A2*p4;
%     p5v = mass_cm3 + A3*p5;
%     p6v = mass_cm3 + A3*p6;
% 
%     pta3 = mass_cm3 + A3*u_ta3;
%     pta2 = mass_cm2 + A2*u_ta2;
% 
%     psol3 = mass_cm3 + A3*u_sol1;
%     psol2 = mass_cm2 + A2*u_sol2;
% 
%    %plot3([p1v(1) p2v(1)],[p1v(2) p2v(2)],[p1v(3) p2v(3)],'r-','LineWidth',2)
%    view(90,0);
%    grid on
%    %axis([-2 2 -1 1 -1 0])
%    axis([-2 2 -0.3 0.3 -1 -0.4])
%    hold on
%    plot3([p3v(1) p4v(1)],[p3v(2) p4v(2)],[p3v(3) p4v(3)],'Color',[0.4 0.6 0.7],'LineWidth',5)
%    plot3([p5v(1) p6v(1)],[p5v(2) p6v(2)],[p5v(3) p6v(3)],'Color',[0.5 0.5 0.5],'LineWidth',5)
% 
%    plot3(pta3(1),pta3(2),pta3(3), 'b*')
%    plot3(pta2(1),pta2(2),pta2(3), 'b*')
%    plot3([pta3(1) pta2(1)],[pta3(2) pta2(2)],[pta3(3) pta2(3)],'b-','LineWidth',1)
% 
%    plot3(psol3(1),psol3(2),psol3(3), 'k*')
%    plot3(psol2(1),psol2(2),psol2(3), 'k*')   
%    plot3([psol3(1) psol2(1)],[psol3(2) psol2(2)],[psol3(3) psol2(3)],'r-','LineWidth',1)
% 
%    hold off
%    pause(0.05)
%    frame = getframe(gcf); %get frame
%    writeVideo(myVideo, frame);
% end    
% close(myVideo)
% %% 
% 
% %plot(t,F_MTC_sol,'red',t,F_MTC_ta,'blue','LineWidth',2);
% plot(t,Tendon_strain_vec*100,'red','LineWidth',2);
% ax=gca;
% ax.FontSize = 18;
% xlim([0 0.2]);
% 
% %xlim([0 600])
% %xticks(0:100:600)
% %ylim([0 1.5])
% %yticks(0:0.25:1.5)
% xlabel('time, s','FontSize',18)
% ylabel('Tendon strain, %','FontSize',18)
% grid on
% 
% %legend({'Soleus','Tibialis anterior',},'FontSize',18)
% 
% figure
% plot(t,q_sol_vec,'red','LineWidth',2);
% % Test if total energy is constant before adding muscle forces
% 
% % Ek = zeros(1,length(t));
% % Ep = zeros(1,length(t));
% % E = zeros(1,length(t));
% % 
% % for kk = 1:length(t)
% %     R31 = y(kk,3);
% %     e01 = y(kk,4);
% %     e11 = y(kk,5);
% %     e21 = y(kk,6);
% %     e31 = y(kk,7);
% %     R32 = y(kk,10);
% %     e02 = y(kk,11);
% %     e12 = y(kk,12);
% %     e22 = y(kk,13);
% %     e32 = y(kk,14);
% %     R33 = y(kk,17);
% %     e03 = y(kk,18);
% %     e13 = y(kk,19);
% %     e23 = y(kk,20);
% %     e33 = y(kk,21);
% %     
% %     Ek(kk) = 0.5*y(kk,22:42)*Mf3(Ix1_,Ix2_,Ix3_,Iy1_,Iy2_,Iy3_,...
% %         Iz1_,Iz2_,Iz3_,e01,e02,e03,e11,e12,e13,e21,e22,e23,e31,...
% %         e32,e33,m1,m2,m3)*y(kk,22:42).';            % 2.85
% %     Ep(kk) = m1*g*(R31) + m2*g*(R32) + m3*g*R33;        % Ep = mgh
% %     E(kk) = Ek(kk) + Ep(kk);
% % end
% % figure(1);
% % plot(t,Ek, t,Ep, t,E)
% % xlabel('Time [s]')
% % ylabel('System energy [J]')
% % legend('Kinetic', 'Potential', 'Total')
% % grid on

function dot_l_MTC_ta = dot_l_mtc_3points(body1coordvec, body2coordvec, dot_body1coordvec, dot_body2coordvec, u1, u2, u3) 

% initial mtc length

% muscle attachment point in body 3
body3_cm = body1coordvec(1:3); 
A3 = Af3(body1coordvec(4),body1coordvec(5),...
    body1coordvec(6),body1coordvec(7));
p1 = body3_cm + A3*u1;

% muscle attachment point in body 2
body2_cm = body2coordvec(1:3);
A2 = Af2(body2coordvec(4),body2coordvec(5),...
    body2coordvec(6),body2coordvec(7));
p2 = body2_cm + A2*u2;

% 3rd attachment point, also in body 2
p3 = body2_cm + A2*u3;


% muscle-tendon length, distance between two points of the bodies
% l_mtc_ta = sqrt((p2(1)-p1(1)).^2 + (p2(2)-p1(2)).^2 + ...
%     (p2(3)-p1(3)).^2) + sqrt((p2(1)-p3(1)).^2 + (p2(2)-p3(2)).^2 + ...
%     (p2(3)-p3(3)).^2);

l_mtc_ta_1 = sqrt((p2(1)-p1(1)).^2 + (p2(2)-p1(2)).^2 + ...
    (p2(3)-p1(3)).^2);
l_mtc_ta_2 = sqrt((p2(1)-p3(1)).^2 + (p2(2)-p3(2)).^2 + ...
    (p2(3)-p3(3)).^2);

l_mtc_ta = l_mtc_ta_1 + l_mtc_ta_2;

% Calculate the Jacobian matrices
J1 = Jac(body1coordvec(4),body1coordvec(5),body1coordvec(6),body1coordvec(7),u1(1),u1(2),u1(3));
J2_1 = Jac(body2coordvec(4),body2coordvec(5),body2coordvec(6),body2coordvec(7),u2(1),u2(2),u2(3));
J2_2 = Jac(body2coordvec(4),body2coordvec(5),body2coordvec(6),body2coordvec(7),u3(1),u3(2),u3(3));

% Calculate the time derivatives of the positions
dp1dt = J1 * dot_body1coordvec;
dp2dt = J2_1 * dot_body2coordvec;
dp3dt = J2_2 * dot_body2coordvec;

% Calculate the time derivative of the distance
dot_l_MTC_ta_1 = dot(dp2dt - dp1dt, (p2 - p1)) / l_mtc_ta_1;
dot_l_MTC_ta_2 = dot(dp2dt - dp3dt, (p2 - p3)) / l_mtc_ta_2;
dot_l_MTC_ta = dot_l_MTC_ta_1 + dot_l_MTC_ta_2;

end
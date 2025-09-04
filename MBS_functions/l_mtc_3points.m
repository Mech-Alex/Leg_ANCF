function l_mtc = l_mtc_3points(body1coordvec, body2coordvec, u1, u2, u3) 

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
l_mtc = sqrt((p2(1)-p1(1)).^2 + (p2(2)-p1(2)).^2 + ...
    (p2(3)-p1(3)).^2) + sqrt((p2(1)-p3(1)).^2 + (p2(2)-p3(2)).^2 + ...
    (p2(3)-p3(3)).^2);

end
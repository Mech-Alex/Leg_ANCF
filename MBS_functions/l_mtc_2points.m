function l_mtc_sol = l_mtc_2points(body1coordvec, body2coordvec, u1, u2) 
% muscle-tendon complex (MTC) length
% calculates length of the vector between attachment points

% SOL attachment point in body 3, heel
body3_cm = body1coordvec(1:3); 
A3 = Af3(body1coordvec(4),body1coordvec(5),...
    body1coordvec(6),body1coordvec(7));
p1_sol = body3_cm + A3*u1;

% SOL attachment point in body 2
body2_cm = body2coordvec(1:3);
A2 = Af2(body2coordvec(4),body2coordvec(5),...
    body2coordvec(6),body2coordvec(7));
p2_sol = body2_cm + A2*u2;

% vector between points
p_sol = p2_sol - p1_sol;
% length of vector
l_mtc_sol = sqrt(sum(p_sol.^2));

end
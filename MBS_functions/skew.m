function Skew_matrix = skew(vector)
%SKEW Returns the skew-symmetric form of the given 3x1 or 1x3 vector
Skew_matrix = [0,-vector(3),vector(2);
               vector(3),0,-vector(1);
               -vector(2),vector(1),0];
end


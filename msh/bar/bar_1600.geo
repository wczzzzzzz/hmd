
a = 8;
n = 1600;

Point(1) = {0, 0, 0};
Point(2) = {a, 0, 0};

Line(1) = {1,2};

Transfinite Curve{1} = n+1;

Physical Point("Γᵗ") = {2};
Physical Point("Γᵍ") = {1};
Physical Curve("Ω") = {1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 1;

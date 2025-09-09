
a = 8;
n = 20;
prog = 1.05;
//c = 0.2;

Point(1) = {0, 0, 0};
Point(2) = {a, 0, 0};
//Point(1) = {0, 0, 0, c};
//Point(2) = {a, 0, 0, c};

Line(1) = {1,2};

//Transfinite Curve{1} = n+1;
Transfinite Curve{1} =  n+1 Using Progression prog;

Physical Point("Γ¹") = {1};
Physical Point("Γ²") = {2};

Physical Curve("Ω") = {1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 1;

a = 1.0;
b = 2.0;
n = 21;
m = 21;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  a, 0.0, 0.0};
Point(3) = {  a,   b, 0.0};
Point(4) = {0.0,   b, 0.0};
Point(5) = {0.0, b/2, 0.0};
Point(6) = {  a, b/2, 0.0};

Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 1};
Line(7) = {5, 6};

Curve Loop(1) = {1, 2, -7, 6};
Curve Loop(2) = {7, 3, 4, 5};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Curve{1, 4, 7} = n;
Transfinite Curve{2, 3, 5, 6} = m;

Transfinite Surface{1}Right;
Transfinite Surface{2};

Physical Curve("Γ¹") = {1};
Physical Curve("Γ²") = {2, 3};
Physical Curve("Γ³") = {4};
Physical Curve("Γ⁴") = {5, 6};
Physical Surface("Ω") = {1, 2};

Mesh.Algorithm = 9;
//Mesh.SecondOrderIncomplete = 1;
Mesh 2;
//SetOrder 2;

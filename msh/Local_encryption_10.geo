
a = 4.0;
b = 4.0;
c1 = 0.05;
c2 = 0.4;

Point(1) = {0.0, 0.0, 0.0, c2};
Point(2) = {  a, 0, 0.0, c2};
Point(3) = {  a,  b, 0.0, c2};
Point(4) = {0.0,  b, 0.0, c2};
Point(5) = {  0,  1.0, 0.0, c1};
Point(6) = {  3.0,  4.0, 0.0, c1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,6};
Line(4) = {6,4};
Line(5) = {4,5};
Line(6) = {5,1};
Line(7) = {5,6};

Curve Loop(1) = {1,2,3,4,5,6};

Plane Surface(1) = {1};

Physical Curve("Γ¹") = {1};
Physical Curve("Γ²") = {2};
Physical Curve("Γ³") = {3,4};
Physical Curve("Γ⁴") = {5,6};
Physical Surface("Ω") = {1};

Curve{7} In Surface{1};
Field[1] = Distance;
Field[1].CurvesList = [5];
Field[1].Sampling = 100;

Mesh.Algorithm = 4;
Mesh.MshFileVersion = 2;
Mesh 2;
//SetOrder 2;
//RecombineMesh;

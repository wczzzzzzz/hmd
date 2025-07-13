
a = 4.0;
b = 8.0;
c = 0.8;
//n = 9;
//m = 16;

Point(1) = {0.0, 0.0, 0.0, c};
Point(2) = {  a, 0.0, 0.0, c};
Point(3) = {  a,   b, 0.0, c};
Point(4) = {0.0,   b, 0.0, c};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

//Transfinite Curve{1,2,3,4} = n;
//Transfinite Curve{2,4} = m;
//Transfinite Curve{1,3} = n;
//Transfinite Curve{2,-4} = 10 Using Progression 1.4;

//Transfinite Surface{1};

Physical Curve("Γ¹") = {1};
Physical Curve("Γ²") = {2};
Physical Curve("Γ³") = {3};
Physical Curve("Γ⁴") = {4};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 2;
//Mesh.Algorithm = 5;
Mesh.MshFileVersion = 2;

Mesh.SecondOrderIncomplete = 1;
Mesh 2;
//RecombineMesh;

//RefineMesh;


//SetOrder 2;

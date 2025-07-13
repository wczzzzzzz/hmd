c1 = 0.002;
c2 = 0.4;

Point(1) = {0, 0, 0, c2};
Point(2) = {2, 0, 0, c2};
Point(3) = {2, 4, 0, c2};
Point(4) = {0, 4, 0, c2};
Point(5) = {0, 0.5, 0, c1};
Point(6) = {1.5, 2, 0, c1};
Point(7) = {2, 2.5, 0, c1};


Line(1) = {1,2};
Line(2) = {2,7};
Line(3) = {7,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,1};
Line(7) = {5,6};
Line(8) = {6,7};




Curve Loop(1) = {5,4,3,8,7};
Curve Loop(2) = {1,2,-8,-7,6};


Plane Surface(1) = {1};
Plane Surface(2) = {2};

//Curve{7,8} In Surface{1};

Physical Curve("Γᵍ₁") = {1};
Physical Curve("Γᵍ₂ ") = {4};
Physical Curve("Γᶜ") = {7};
Physical Surface("Ω") = {1,2};


Mesh.Algorithm = 1;
Mesh.MshFileVersion = 3;
Mesh 2;
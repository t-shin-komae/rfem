SetFactory("OpenCASCADE");
lc = 0.04;
//+
Point(1) = {0,0,0,lc};
//+
Point(2) = {1,0,0,lc};
//+
Point(3) = {1,1,0,lc};
//+
Point(4) = {0,1,0,lc};
//+
Point(6) = {0.7, 0.5, 0, lc/2};
//+
Point(7) = {0.5, 0.7, -0, lc/2};
//+
Point(8) = {0.3, 0.5, -0, lc/2};
//+
Point(9) = {0.5, 0.3, 0, lc/2};
//+
Point(10) = {0.3, 0.4, -0, lc/2};
//+
Point(11) = {0.3, 0.6, -0, lc/2};
//+
Point(12) = {0.7, 0.6, 0, lc/2};
//+
Point(13) = {0.7, 0.4, 0, lc/2};
//+
Point(14) = {0.7, 0.3, 0, lc/2};
//+
Point(15) = {0.3, 0.3, -0, lc/2};
//+
Point(16) = {0.3, 0.7, -0, lc/2};
//+
Point(17) = {0.7, 0.7, 0, lc/2};
//+
Point(18) = {0.6, 0.7, 0, lc/2};
//+
Point(19) = {0.4, 0.7, -0, lc/2};
//+
Point(20) = {0.4, 0.3, -0, lc/2};
//+
Point(21) = {0.6, 0.3, 0, lc/2};
//+
Line(1) = {1,2};
//+
Line(2) = {2,3};
//+
Line(3) = {3,4};
//+
Line(4) = {4,1};
//+
Circle(5) = {0.5, 0.5, 0, 0.2, 0, 2*Pi};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("IN") = {1};
//+
Physical Curve("OUT") = {3};
//+
Physical Curve("WALL") = {4, 5, 2};
//+
Physical Surface("Domain") = {1};

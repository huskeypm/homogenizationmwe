lc = 0.1;  
lco = 1.;
del = 5.;
eps = 0.1;

// outer 
Point(1) = {10, 10, -0, lco};
Point(2) = {0, 10.0, -0, lco};
Point(3) = {0, 0, -0, lco};
Point(4) = {10, 0, -0, lco};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(1) = {3, 4, 1, 2};


/// Inclusions 
// inner
Point(5) = {1, 1, 0,lc};
Point(6) = {1, 1+del,0, lc};
Point(7) = {1+del, 1+del, 0, lc};
Point(8) = {1+del, 1,    0, lc};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
// outer
Point(9) = {1+eps, 1+eps, 0,lc};
Point(10) = {1+eps, 1+del-eps,0, lc};
Point(11) = {1+del-eps, 1+del-eps, 0, lc};
Point(12) = {1+del-eps, 1+eps,    0, lc};
Line(9) = {12, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line Loop(13) = {6, 7, 8, 5};
Plane Surface(14) = {1, 13};
Line Loop(15) = {11, 12, 9, 10};
Plane Surface(16) = {15};

lc=0.02;
Point(5) = {0.0, 0.0, -1, lc };
Point(6) = {0.0,-0.05,-1 , lc };
Point(7) = {0.0, 0.05,-1, lc };
Point(8) ={-0.05,0.0, -1, lc };
Point(9) = {0.05, 0.0,-1, lc };
Circle(1) = {8, 5, 7};
Circle(2) = {7, 5, 9};
Circle(3) = {9, 5, 6};
Circle(4) = {6, 5, 8};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
Extrude {0, 0, 2} {
  Surface{6};
}

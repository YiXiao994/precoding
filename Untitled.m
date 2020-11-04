gpvar x y z;

m = x + y + 2*z;

constraints = [m <=x, z <= y, y/(x+z) >= 5];
[obj_value, solution, status] = gpsolve(m,constraints);
assign(solution);
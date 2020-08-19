function X = inverse_map_Q4(coords,x)

X = [0;0];
N = [0.25*(1+X(1))*(1+X(2));
    0.25*(1-X(1))*(1+X(2));
    0.25*(1-X(1))*(1-X(2));
    0.25*(1+X(1))*(1-X(2))];
xn = coords*N;
delta = x-xn;
converged = (norm(delta) < 1e-12);

itr = 0;

while ~converged
    itr = itr+1;
    
    dN_dxi = [0.25*(1+X(2));
        -0.25*(1+X(2));
        -0.25*(1-X(2));
        0.25*(1-X(2))];
    dN_deta = ...
        [0.25*(1+X(1));
        0.25*(1-X(1));
        -0.25*(1-X(1));
        -0.25*(1+X(1))];
    J = [coords(1,:)*dN_dxi,coords(1,:)*dN_deta;
        coords(2,:)*dN_dxi,coords(2,:)*dN_deta];
    delta_X = J\delta;
    X = X+delta_X;
    
    N = [0.25*(1+X(1))*(1+X(2));
        0.25*(1-X(1))*(1+X(2));
        0.25*(1-X(1))*(1-X(2));
        0.25*(1+X(1))*(1-X(2))];
    xn = coords*N;
    delta = x-xn;
    converged = (norm(delta) < 1e-12);
    
    fprintf('iteration %d, residual = %.3E\n',itr,norm(delta));
end
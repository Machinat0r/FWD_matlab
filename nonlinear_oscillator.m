function dxdt = nonlinear_oscillator(t, x)
    omega = 2*pi;
    dxdt = zeros(2,1);
    dxdt(1) = x(2);
    dxdt(2) = -sin(x(1)) - 0.1*x(2) + 0.5*cos(omega*t);
end

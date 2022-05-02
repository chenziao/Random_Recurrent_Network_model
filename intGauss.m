function I = intGauss( f, r, w )
% Integrate function f with respect to standard Gaussian variable
% I = 1/sqrt(2*pi) * int f(x)e^(-x^2/2)
% r - roots of Hermite polynomial
% w - associated weigths
I = sum(f(2^0.5*r).*w)/pi^0.5;
end
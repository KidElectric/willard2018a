function y=logmodulus(x)
% y=logmodulus(x)
% perform log modulus of x: y=sign(x) .* log10(abs(x)+1)

y=sign(x) .* log10(abs(x)+1);% 
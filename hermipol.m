function [p,r,w] = hermipol(n)
p = {1,[2 0]};
for k = 2:n
   p{k+1} = [2*p{k},0] - 2*(k-1)*[0,0,p{k-1}];
end
if nargout>1
    r = sort(real(roots(p{n+1})))';
end
if nargout>2
	w = (2^(n-1)/n^2*pi^0.5)*(factorial(n)./polyval(p{n},r).^2);
end
p = p{n+1};
end

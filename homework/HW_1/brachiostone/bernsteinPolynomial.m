function B = bernsteinPolynomial(p,domain,x)
scaleFactor = 1/((domain(2)-domain(1))^p);
B = zeros(length(x),p+1);
for a=1:p+1
    binomialCoefficient = factorial(p)/(factorial(a-1)*factorial(p+1-a));
    B(:,a) = scaleFactor * binomialCoefficient * (domain(2)-x).^(p-(a-1)) .*(-domain(1)+x).^(a-1);
end
end
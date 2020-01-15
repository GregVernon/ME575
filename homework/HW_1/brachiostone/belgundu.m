function time = belgundu(y,dx,m,g,mu)

time = zeros(length(y)-1,1);
v0 = 0;
for ii = 1:length(y)-1
    if ii == 1
        y0 = 1;
    else
        y0 = y(ii);
    end
    
    if ii == length(y)-1
        y1 = 0;
    else
        y1 = y(ii+1);
    end
    
    theta = atan2((y0-y1),dx);
    
    L = sqrt(dx^2 + (y1-y0)^2);
    v1 = sqrt((2/m) * ((-mu * g * cos(theta) * L) - (m * g * (y1 - y0))) + v0^2);
    time(ii) = ((v1 - v0)*L / (g*(y0-y1)));
    v0 = v1;
end

time = sum(time);
end

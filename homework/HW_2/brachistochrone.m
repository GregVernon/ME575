function [f,g] = brachistochrone(yint)

    global fcalls
    if fcalls > 1e4
        return
    end

    mu_k = 0.3;

    y = cat(1,1, yint, 0);
    n = length(y);
    x = linspace(0.0, 1.0, n);
    g = zeros([n-2,1]);

    T = 0.0;
    for i = 1:n-1 
        ds = sqrt((x(i+1) - x(i))^2 + (y(i+1) - y(i))^2);

        if 1 - y(i+1) - mu_k*x(i+1) < 0 || 1 - y(i) - mu_k*x(i) < 0
            T = T + 10;

        else

            vbar = sqrt(1 - y(i+1) - mu_k*x(i+1)) + sqrt(1 - y(i) - mu_k*x(i));

            % gradient
            if i > 1
                dsdyi = 0.5/ds*2*(y(i+1) - y(i)) * -1;
                dvdyi = 0.5/sqrt(1 - y(i) - mu_k*x(i)) * -1;
                dtdyi = (vbar*dsdyi - ds*dvdyi)/(vbar^2);
                g(i-1) = g(i-1) + dtdyi;
            end
            if i < n-1
                dsdyip = 0.5/ds*2*(y(i+1) - y(i));
                dvdyip = 0.5/sqrt(1 - y(i+1) - mu_k*x(i+1)) * -1;
                dtdyip = (vbar*dsdyip - ds*dvdyip)/(vbar^2);
                g(i) = g(i) + dtdyip;
            end
            
            T = T + ds/vbar;
        end
    end
    
    f = T;
    
    fcalls = fcalls + 1;

end          


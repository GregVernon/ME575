clear

r1 = linspace(0,1,100);
r2 = linspace(0,1,100);
theta = linspace(0,pi,100);
[r1,r2,theta] = meshgrid(r1,r2,theta);
f = zeros(size(r1));

for n = 1:numel(r1)
    f(n) = computeFreq(r1(n), r2(n), theta(n));
end

x1 = r1;
y1 = zeros(size(r1));

x2 = r2.*cos(theta);
y2 = r2.*sin(theta);

dR = sqrt((x2-x1).^2 + (y2-y1).^2);
rMin = min(r1,r2);
rMax = max(r1, r2);

fid = fopen("data_3.csv",'w+');
fprintf(fid, "r1\tr2\ttheta\tfreq\tdR\trMin\trMax\n");
fprintf(fid,"%12.8f \t %12.8f \t %12.8f \t %12.8f \t %12.8f \t %12.8f \t %12.8f \n",[r1(:) r2(:) theta(:) f(:) dR(:) rMin(:) rMax(:)]');
fclose(fid);



clear
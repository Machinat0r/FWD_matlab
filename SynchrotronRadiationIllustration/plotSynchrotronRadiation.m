figure('Color','white');

% plot blue vectors to be the magnetic field
for zm = -9:9:9
    for ym = -9:9:9
        z = [-15:1:70]';
        npts = length(z);
        x = ym.*ones(npts,1);
        y = zm.*ones(npts,1);
        width = ones(npts,1);
        n = 10;
        verts = [x y z];
        [vx,vy,vz]=spaghetti_XYZ(verts,width,n);
        surface(vx,vy,vz,'facecolor','b','edgecolor','none');
        hold on;
        arrow3D([x(end) y(end) z(end)], [0 0 3.5], 'b', 0.1, 0.5);
    end
end

% helix/spiral path of the charged particle
r = 20;
a = 0.04;
theta = [0:360*4]';
hx = r.*cosd(theta);
hy = r.*sind(theta);
hz = a.*theta;
np = length(theta);
colorCode = [1 0.5 0];
C = [ones(np,1)' zeros(np,1)' zeros(np,1)'];
width = 2.0*ones(np,1);
n = 20;
verts = [hx hy hz];
[vx,vy,vz]=spaghetti_XYZ(verts,width,n);
s = surf(vx,vy,vz, 'FaceColor', colorCode, 'EdgeColor', 'none');
hold on;
% cone/arrow end to the spiral showing its direction of motion
theta2 = [360*4:360*4+20]';
hx2 = r.*cosd(theta2);
hy2 = r.*sind(theta2);
hz2 = a.*theta2;
np2 = length(theta2);
C2 = [ones(np2,1)' zeros(np2,1)' zeros(np2,1)'];
width2 = 4:-0.2:0;
n2 = 20;
verts2 = [hx2 hy2 hz2];
[vx2,vy2,vz2]=spaghetti_XYZ(verts2,width2,n);
s2 = surf(vx2,vy2,vz2, 'FaceColor', colorCode, 'EdgeColor', 'none');
hold on;
% yellow sphere to represent the charged particle
[spx,spy,spz] = sphere; 
spr = 3;
surf(spx.*spr+hx(1),spy.*spr+hy(1),spz.*spr+hz(1), ...
    'FaceColor', colorCode, 'EdgeColor', 'none');

% photon that is emitted
t = 1:360*5;
phx = (-24.*t./3600.0)+8;
phy = (hy(802)+3).*ones(1,length(t)) + 15.*t./3600.0; 
phz = (+2.*sind(t)) + 4.*t./3600.0 + hz(802);
verts = [phx' phy' phz'];
width = 0.5.*ones(length(t),1);
[vx,vy,vz]=spaghetti_XYZ(verts,width,n);
s = surf(vx,vy,vz, 'FaceColor', 'r', 'EdgeColor', 'none');
arrow3D([phx(end) phy(end) phz(end)], [-5 3 1], 'r', 0.2, 0.5);

% pretty shading and perspective
view(3);
view(104,22);
axis off;
xlabel('x');
ylabel('y');
grid off;
axis equal;
light;
lighting phong;
camlight('left');

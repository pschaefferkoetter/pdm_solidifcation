function plot_mpf (f1, coord, phi, step)

if (nargin == 3)
    step = 0;
end

figure(f1)
tri = delaunay(coord(1,:), coord(2,:));
trisurf(tri, coord(1,:), coord(2,:), max(phi(:,2:end),[],2)');
set(gca,'dataaspectratio', [1 1 500])
title(['Step ',num2str(step)])
view(0,90)
shading interp
caxis([0 1])
colormap jet
colorbar
drawnow

end
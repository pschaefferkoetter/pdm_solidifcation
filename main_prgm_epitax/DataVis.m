function [] = DataVis(figure_num,pos,var,set_c)

[m, n] = size(pos);

if m < n
    pos = pos';
end
    
x = pos(:,1) ;
y = pos(:,2) ;
z = var;
% Grid 
x0 = min(x) ; x1 = max(x) ;
y0 = min(y) ; y1 = max(y) ;
N = 100;
xl = linspace(x0,x1,N) ; 
yl = linspace(y0,y1,N) ; 
[X,Y] = meshgrid(xl,yl) ;
%do inteprolation 
P = [x,y] ; V = z ;
F = scatteredInterpolant(P,V) ;
F.Method = 'natural';
F.ExtrapolationMethod = 'none' ;  % none if you dont want to extrapolate
% Take points lying insuide the region
pq = [X(:),Y(:)] ; 
vq = F(pq) ;
Z = vq ;
Z = reshape(Z,size(X)) ;

figure(figure_num)
surf(X,Y,Z,'edgecolor','interp');
colorbar
colormap jet
shading interp
view([0 90])


if set_c~='n'
    caxis([set_c(1),set_c(2)])
else
    caxis('auto')
end
end


clear; close all; clc;

%Define domain limits
L =1.0;

%Define the dimension
dim = 2;

%Define minimum spacing increment
h_min = L/20;

%Define maximum spacing increment
h_max = L/10;

%Define munber of points within the domain
n_pts = 2500;

%Preallocate memory and intitlize 
coord = zeros(dim, n_pts);

%Define Start Point
coord(:,1) = L/2;
x_old = coord(:,1);

%Define the domain corners
coord(:,2) = [0;0];
coord(:,3) = [L;0];
coord(:,4) = [L;L];
coord(:,5) = [0;L];

%Start Counter 
i_pnt = 6;

%Take a walk around the park
while i_pnt <= n_pts
    
    %Generate a random direction
    m = cos(rand*2*pi);

    Direction = [m; sqrt(1 - m^2)]; 
    
    %Generate a random legth within limits 
    S = h_min + rand*(h_max - h_min);
    
    %Update position
    x_new = x_old + S*Direction;
    
    %Is x_new too close to another point or outside of the domain:
    dist2pnts = vecnorm(x_new - coord(:,1:i_pnt-1));
    mindist2pnts = min(dist2pnts);
    
    if x_new(1) > L || x_new(2) > L || x_new(1) < 0 || x_new(2) < 0 || mindist2pnts < h_min
       continue
    else 
        coord(:,i_pnt) = x_new;
        x_old = x_new;
        
        figure(1)
        scatter(coord(1,1:i_pnt),coord(2,1:i_pnt))
        pause(.01)
        
        i_pnt = i_pnt + 1;
       
    end 
    
    
end    
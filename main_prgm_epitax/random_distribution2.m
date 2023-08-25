clear; close all; clc;

%Define domain limits
L =1.0;

%Define the dimension
dim = 2;

%Define minimum spacing increment
h_min = .008;

%Define maximum spacing increment
h_max = L/10;

%Define munber of points within the domain
n_pts = 10000;

%Preallocate memory and intitlize 
coord = zeros(dim, n_pts);


%Define the domain corners
coord(:,1) = [0;0];
coord(:,2) = [L;0];
coord(:,3) = [L;L];
coord(:,4) = [0;L];

current_pnt = 5;
no_pnt_cntr = 0;

tic

while current_pnt <= n_pts
    

    previous_pnt = current_pnt - 1;
    current_loc = [rand; rand];
    
    mindist2pnts = min(vecnorm(current_loc - coord(:,1:previous_pnt)));
    
    if mindist2pnts < h_min
        no_pnt_cntr = no_pnt_cntr + 1;
        if no_pnt_cntr == 10000
            disp('Thats the best your going to get!')
            previous_pnt
            break
        end 
        continue
    else
        coord(:,current_pnt) = current_loc;
        current_pnt = current_pnt + 1;
        no_pnt_cntr = 0;
    end 
      
%     figure(1)
%     scatter(coord(1,1:current_pnt),coord(2,1:current_pnt))
%     pause(.01)
end
toc
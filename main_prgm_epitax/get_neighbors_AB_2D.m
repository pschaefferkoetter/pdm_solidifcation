function [nbr_info, coord, n_node] = get_neighbors_AB_2D (coord, Length, min_dist, perbc, nnbr, TEST)
% -------------------------------------------------------------------------
% Compute nbr info for calculation of PDM differential operators
% -------------------------------------------------------------------------
% INPUTS:
% - coord (2 x n_node) list of coordinates
% - box (1 x 4 or 4 x 1) vector containing [xmin xmax ymin ymax] of box.
%       If domain is not a rectangle, perbc should be off and box can be
%       set to anything desired because it will not enter into the
%       calculations in this function.
% - min_dist (1 x 1) minimum distance between nodes.
% - perbc (1 x 2 or 2 x 1) vector identifying whether perbc are to be used
%       in each direction (0:no perbc, 1:perbc).  e.g. [1 0] means perbc in
%       x but not y, [1 1] means perbc in both directions
% - nnbr (1 x 1) constant number of neighbors specified for knnsearch
% - TEST (1 x 1) optional argument telling function whether to run
%       graphical test on node neighbors. 0 (off) by default.
% -------------------------------------------------------------------------
% OUTPUTS:
% - nbr_info (struct).  Members:
%   - nblist (n_node x 1 cell of 1 x nnbr vectors)
%       neighbor ID lists for each node
%   - dilp (n_node x 1 cell of double values)
%       dilation parameter (distance to farthest neighbor in compact
%       support) for each node
%   - dist_list (n_node x 1 cell of 1 x nnbr vectors)
%       lists of distances between each node and each of its compact
%       support neighbors. Adjustment for periodicity built in.
%   - xixj_list (n_node x 1 cell of 2 x nnbr matrices)
%       lists of vectors (xj-xi) between each node and each of its compact
%       support neighbors for use in get_pvec2d. Adjustment for periodicity
%       built in.
% - coord (2 x n_node) updated list of coordinates with top and/or right
%       edge (or neither) removed depending on perbc
% - n_node (1 x 1) new number of nodes after removing nodes
% -------------------------------------------------------------------------

if (nargin == 5)
    TEST = 0;
end

disp('Starting get_neighbors')
start_time = tic;

% Get total number of nodes and box dimensions
n_node = length(coord);
Lx = box(2) - box(1);
Ly = box(4) - box(3);

% Throw out coordinates on right edge, if necessary
TOL = 0.01*min_dist;
keep = ones(1,length(coord));
if (perbc(1))
    for i = 1:n_node
        if ( abs(coord(1,i) - box(2)) < TOL )
            keep(i) = 0;
        end
    end
end
% Throw out coordinates on top edge, if necessary
if (perbc(2))
    for i = 1:n_node
        if ( abs(coord(2,i) - box(4)) < TOL )
            keep(i) = 0;
        end
    end
end
coord = coord(:,logical(keep));
n_node = length(coord);

% Make clones of coord as needed
if ( perbc(1) && perbc(2) )
    ctr = 1;
    coord_expand = zeros(2,9*n_node);
    % Original
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:);
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:);
    ctr = ctr + 1;
    % xmin, ymin
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:) - Lx;
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:) - Ly;
    ctr = ctr + 1;
    % x0, ymin
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:);
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:) - Ly;
    ctr = ctr + 1;
    % xmax, ymin
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:) + Lx;
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:) - Ly;
    ctr = ctr + 1;
    % xmax, y0
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:) + Lx;
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:);
    ctr = ctr + 1;
    % xmax, ymax
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:) + Lx;
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:) + Ly;
    ctr = ctr + 1;
    % x0, ymax
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:);
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:) + Ly;
    ctr = ctr + 1;
    % xmin, ymax
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:) - Lx;
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:) + Ly;
    ctr = ctr + 1;
    % xmin, y0
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:) - Lx;
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:);
elseif (  perbc(1) && ~perbc(2) )
    ctr = 1;
    coord_expand = zeros(2,3*n_node);
    % Original
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:);
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:);
    ctr = ctr + 1;
    % xmin, y0
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:) - Lx;
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:);
    ctr = ctr + 1;
    % xmax, y0
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:) + Lx;
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:);
elseif ( ~perbc(1) &&  perbc(2) )
    ctr = 1;
    coord_expand = zeros(2,3*n_node);
    % Original
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:);
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:);
    ctr = ctr + 1;
    % x0, ymin
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:);
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:) - Ly;
    ctr = ctr + 1;
    % x0, ymax
    coord_expand(1,1+(ctr-1)*n_node:ctr*n_node) = coord(1,:);
    coord_expand(2,1+(ctr-1)*n_node:ctr*n_node) = coord(2,:) + Ly;
else
    % perbc is off for both x and y
    coord_expand = coord;
end

% Build kd search tree
Mdl = KDTreeSearcher(coord_expand');
nblist = cell(n_node,1);
dist_list = cell(n_node,1);
% Search for k nearest neighbors
for i = 1:n_node
    [nblist{i}, dist_list{i}] = knnsearch(Mdl,coord(:,i)','k',nnbr);
end

% Compute vectors between nodes
% Modify neighbor ids using mod() 
% Get dilp
xixj_list = cell(n_node,1);
dilp = zeros(n_node,1);
for i = 1:n_node
    xixj_list{i} = zeros(2,nnbr);
    for j = 1:length(nblist{i})
        jj = nblist{i}(j);
        xixj_list{i}(:,j) = coord_expand(:,jj) - coord(:,i);
        nblist{i}(j) = 1+mod(jj-1,n_node);
        if jj > n_node
            [jj nblist{i}(j)]
        end
    end
    dilp(i) = dist_list{i}(end);
end

% Place variables in output data structure
nbr_info = struct;
nbr_info.nblist = nblist;
nbr_info.dilp = dilp;
nbr_info.dist_list = dist_list;
nbr_info.xixj_list = xixj_list;

disp(['get_neighbors finished in ',num2str(toc(start_time)),' sec']);

if (TEST)
    figure
    for i = 1:n_node
        % Plot all coord
        plot(coord(1,:), coord(2,:), 'k.')
        set(gca, 'dataaspectratio', [1 1 1])
        hold on
        % Plot neighbors of coord(i)
        plot(coord(1,nblist{i}), coord(2,nblist{i}), 'bo')
        hold on
        % Plot vectors from node i to each neighbor
        quiver(coord(1,i) + zeros(1,nnbr), coord(2,i) + zeros(1,nnbr), xixj_list{i}(1,:), xixj_list{i}(2,:), 'k')
        xlim([box(1)-0.1, box(2)+0.1])
        ylim([box(3)-0.1, box(4)+0.1])
        hold on
        % Plot circle around center node with radius dilp(i)
        plot(coord(1,i) + dilp(i)*cos(0:pi/100:2*pi), coord(2,i) + dilp(i)*sin(0:pi/100:2*pi), 'r', 'LineWidth', 1.0)
        pause(0.5)
        hold off
    end
end

end
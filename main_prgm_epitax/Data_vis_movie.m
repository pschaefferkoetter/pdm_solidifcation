function [void] = Data_vis_movie(STORE_var,fld_num,coord,step_range)

if step_range ~= 'n'
    start = step_range(1);
    stop  = step_range(2);
else 
    start = 1;
    stop  = n_step;
end


[n_step,~] = size(STORE_var);

switch fld_num
    case 2
        %cmin = min(STORE_var{start,fld_num}(:,2));
        cmin = 0.000;
        cmax = max(STORE_var{stop,fld_num}(:,2));
    otherwise
        %cmin = min(STORE_var{start,fld_num}(:,1));
        cmin = 0.95;
        %cmax = max(STORE_var{stop,fld_num}(:,1));
        cmax = 1.0;
end
    

for i_step = start:stop
    switch fld_num
        case 2
            fld = max(STORE_var{i_step,fld_num}(:,2:end-1),[],2);
        otherwise
            fld = 1-STORE_var{i_step,fld_num};
    end
    figure(1)
    DataVis(1,coord,fld,[cmin,cmax])
    set(gca,'XColor', 'none','YColor','none')
    set(gcf,'color','w')
    pause( 0.001 );
    
end

void = 0;
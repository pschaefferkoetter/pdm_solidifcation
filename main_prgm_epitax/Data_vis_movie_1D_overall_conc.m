function [void] = Data_vis_movie_1D_overall_conc(STORE,coord,plt_nodes,step_range)

[n_step,~] = size(STORE);
c_max = max(STORE{250,5}(plt_nodes));

if step_range ~= 'n'
    start = step_range(1);
    stop  = step_range(2);
else 
    start = 1;
    stop  = n_step;
end

for i_step = start:stop

    figure(1)
    subplot(2,1,1)
    plot(coord(1,plt_nodes),STORE{i_step,5}(plt_nodes),'b')
    ylim([0,c_max])  
%      xticks(coord(1,plt_nodes))
%      xticklabels(coord(1,plt_nodes))
    title('Overall Concentration')

    subplot(2,1,2)
    plot(coord(1,plt_nodes),max(STORE{i_step,2}(plt_nodes,2:end),[],2),'r')
    ylim([0,1])  
    xticks(coord(1,plt_nodes))
% %     xticklabels(coord(1,plt_nodes))
% %     ylim([0,1])
    title('Solid Phase Order Parameter')
    
    
    pause( 0.001 );    
    
    
end

void = 0;
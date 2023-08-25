function [void] = plot_cross_section(coord,var_field,phase_indx,plt_nodes)




    plot(coord(1,plt_nodes),var_field(plt_nodes,phase_indx))
    xticks(coord(1,plt_nodes))
    xticklabels(plt_nodes)
    title(['Crossection of solution field for phase ', num2str(phase_indx)])

void = 0;

end
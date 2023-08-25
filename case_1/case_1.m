% Description
%   Runs corresponding case_1_in.m input filed and simulates solidification 
%   for a computing cluster using the corresponding case1.sh file
clear all; clc; close all;

run([mfilename,'_in.m'])
addpath('../../main_prgm_epitax')

run main1_epitax

save('my_data');

%Create ouptut folder
dir_out = ['../../output_folder/',mfilename,'_out'];
mkdir(dir_out);

%Save workspace variables into the folder
save([dir_out,'/STORE'],'STORE','-v7.3');
%save([dir_out,'/plt_nodes'],'plt_nodes');
save([dir_out,'/coord'],'coord');

pnt_txt  = sprintf('%d',num_pts);
grn_txt  = sprintf('%d',num_seed_terms);
poly_txt = sprintf('%d',poly_order);
tstep_txt = sprintf('%0.2e',delta_t);
pf_crit_hgt_txt = sprintf('%0.2f',pf_crit_hight);
num_seed_txt  = sprintf('%d',num_seed_terms);

%Generate text file with relavent data

fid = fopen([dir_out,'/results.txt'],'wt');
fprintf(fid,[mfilename,'\n']);
fprintf(fid,['Discretization Type = ',discretization_type, '\n']);
fprintf(fid,['Number of Points = ',pnt_txt,'\n']);
fprintf(fid,['Number of Grains = ',grn_txt,'\n']);
fprintf(fid,['Order of Polynomial = ',poly_txt,'\n']);
fprintf(fid,['Time Step  = ',tstep_txt,'\n']);
fprintf(fid,['PF Threshold Height  = ',pf_crit_hgt_txt,'\n']);
fprintf(fid,['Number of Seeds = ',num_seed_txt,'\n']);
fclose(fid);
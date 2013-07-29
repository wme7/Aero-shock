% stuDG_driver.m

% script to prepare input file, invoke DG executable, read the output file
% and process results for DG code

clear
clc
close('all')

nel = 2;
nop = 2;
exact_integration = 'F';
ntime = 1000;
time_final = 1.0;
icase = 2;
kstages = 3;

% --- delete old input and output files -----
s = system('rm *.txt');

% --- open input file  ---------------------------
input_file_name = 'dg_input.txt';
infile_fid=fopen(input_file_name,'w');

% ---- write input file data ---------------------
fprintf(infile_fid,'%d \n',nel);
fprintf(infile_fid,'%d \n',nop);
fprintf(infile_fid,'%s \n',exact_integration);
fprintf(infile_fid,'%d \n',ntime);
fprintf(infile_fid,'%g \n',time_final);
fprintf(infile_fid,'%d \n',icase);
fprintf(infile_fid,'%d \n',kstages);

% ---- close the input file when done ---------------

fclose(infile_fid);

% ----- invoke the DG Program ----------------------
tic;
s = system('./stuDG_F90');

program_run_time = toc;

q0 = load('dg_output.txt'); 

% ---- visualization/analysis code ---------------------------
function plot_iter

%Dirichlet case
%resid_1=[];

% periodic case, using Jacobi omega=0.5 to control
% num_MG_Vcycles=1
resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Optimal_Vcycle_Per_it/2d_periodic_1VcylePerIt.dat');
% num_MG_Vcycles=2
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Optimal_Vcycle_Per_it/2d_periodic_2VcylePerIt.dat');
% num_MG_Vcycles=4
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Optimal_Vcycle_Per_it/2d_periodic_4VcylePerIt.dat');
% num_MG_Vcycles=10
resid_4=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Optimal_Vcycle_Per_it/2d_periodic_10VcylePerIt.dat');


% periodic case, based on Gauss-Seidel iteratoin.
% num_MG_Vcycles=1
resid_5=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Optimal_Vcycle_Per_it/2d_dirichlet_1VcylePerIt.dat');
% num_MG_Vcycles=2
resid_6=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Optimal_Vcycle_Per_it/2d_dirichlet_2VcylePerIt.dat');
% num_MG_Vcycles=3
resid_7=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Optimal_Vcycle_Per_it/2d_dirichlet_3VcylePerIt.dat');
% num_MG_Vcycles=4
resid_8=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Optimal_Vcycle_Per_it/2d_dirichlet_4VcylePerIt.dat');



figure 

% peiriodic case
plot([0;resid_5(:,1)], [0;log10(resid_5(:,2))],  'b*');
hold on;
plot([0;resid_6(:,1)], [0;log10(resid_6(:,2))],  'r*'); 
hold on;
plot([0;resid_7(:,1)], [0;log10(resid_7(:,2))],  'k*');
hold on;
plot([0;resid_8(:,1)], [0;log10(resid_8(:,2))],  'c*'); 
hold on;


 plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b^'); 
 hold on;
 plot([0;resid_2(:,1)], [0;log10(resid_2(:,2))],  'r^');
 hold on;
 plot([0;resid_3(:,1)], [0;log10(resid_3(:,2))],  'k^');
 hold on;
 plot([0;resid_4(:,1)], [0;log10(resid_4(:,2))],  'c^');
 hold on;

axis square 
axis([0 70 -12.8 0]);

legend('Diri 1 MG/it', 'Diri 2 MG/it', 'Diri 3 MG/it', 'Diri 4 MG/it', 'Peri 1 MG/it', 'Peri 2 MG/it', 'Peri 4 MG/it', 'Peri 10 MG/it', 'Location','NorthEast')
%text(20, -12.98, 'total no. of MG Sweeps');
%legend('Peri 1 MG/it', 'Peri 2 MG/it', 'Peri 4 MG/it', 'Peri 10 MG/it', 'Location','NorthEast')
%text(20, -13.98, 'total no. of MG Sweeps');


hold on;
text(-12, -5, 'log(r_n/r_0)');


% 3D case. 
% Dirichlet case, not done yet.


disp('3D case');

% periodic case, using Jacobi omega=0.5 to control
% num_MG_Vcycles=1
resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/3d_periodic_1VcylePerIt.dat');
% num_MG_Vcycles=2
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/3d_periodic_2VcylePerIt.dat');
% num_MG_Vcycles=4
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/3d_periodic_4VcylePerIt.dat');
% num_MG_Vcycles=10
resid_4=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/3d_periodic_10VcylePerIt.dat');


% periodic case, based on Gauss-Seidel iteratoin.
% num_MG_Vcycles=1
resid_5=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/3d_dirichlet_1VcylePerIt.dat');
% num_MG_Vcycles=2
resid_6=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/3d_dirichlet_2VcylePerIt.dat');
% num_MG_Vcycles=3
resid_7=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/3d_dirichlet_3VcylePerIt.dat');
% num_MG_Vcycles=4
resid_8=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/3d_dirichlet_4VcylePerIt.dat');

    

figure
%Dirichlet case
plot([0;resid_5(:,1)], [0;log10(resid_5(:,2))],  'b*'); 
hold on;
plot([0;resid_6(:,1)], [0;log10(resid_6(:,2))],  'r*');
hold on;
plot([0;resid_7(:,1)], [0;log10(resid_7(:,2))],  'k*');
hold on;
plot([0;resid_8(:,1)], [0;log10(resid_8(:,2))],  'c*');
hold on;

% periodic case 1-4 Jacobi iterataion, 5-8, GS iteration
 plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b^'); 
 hold on;
 plot([0;resid_2(:,1)], [0;log10(resid_2(:,2))],  'r^');
 hold on;
 plot([0;resid_3(:,1)], [0;log10(resid_3(:,2))],  'k^');
 hold on;
 plot([0;resid_4(:,1)], [0;log10(resid_4(:,2))],  'c^');
 hold on;


axis square 

%legend('Diri 1 MG/it', 'Diri 2 MG/it', 'Diri 3 MG/it', 'Diri 4 MG/it', 'Peri 1 MG/it', 'Peri 2 MG/it', 'Peri 4 MG/it', 'Peri 10 MG/it', 'Location','NorthEast')
%text(20, -12.98, 'total no. of MG Sweeps');

legend('Peri 1 MG/it', 'Peri 2 MG/it', 'Peri 4 MG/it', 'Peri 10 MG/it', 'Location','NorthEast')
text(20, -12.98, 'total no. of MG Sweeps');

hold on;
text(-12, -5, 'log(r_n/r_0)');
axis([0 80 -11.8 0]);



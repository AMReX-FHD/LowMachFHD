function compare_P2P3_sign
figure  % Dirichlet using different V-Cycle numbers

 % all Dirichlet, use 1Vcyle/it 
resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_1VcylePerIt_finegrid_prec2.dat');  %  
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_1VcylePerIt_finegrid_prec3.dat'); %
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_1VcylePerIt_finegrid_prec2_positiveSchur.dat'); %
resid_4=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_1VcylePerIt_finegrid_prec3_positiveSchur.dat'); %


 % all Dirichlet, use 25Vcyle/it 
resid_5=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_25VcylePerIt_finegrid_prec2.dat');  %  
resid_6=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_25VcylePerIt_finegrid_prec3.dat'); %
resid_7=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_25VcylePerIt_finegrid_prec2_positiveSchur.dat'); %
resid_8=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_25VcylePerIt_finegrid_prec3_positiveSchur.dat'); %

%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
plot([0;resid_1(1:end,1)/1], [0;log10(resid_1(1:end,2))],  'b*'); 
hold on;
plot([0;resid_2(1:end,1)/1], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)/1], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;
plot([0;resid_4(1:end,1)/1], [0;log10(resid_4(1:end,2))],  'c*');  
hold on;

plot([0;resid_5(1:end,1)/25], [0;log10(resid_5(1:end,2))],  'bs'); 
hold on;
plot([0;resid_6(1:end,1)/25], [0;log10(resid_6(1:end,2))],  'rs'); 
hold on;
plot([0;resid_7(1:end,1)/25], [0;log10(resid_7(1:end,2))],  'ks');  
hold on;
plot([0;resid_8(1:end,1)/25], [0;log10(resid_8(1:end,2))],  'cs');  
hold on;

%hold on;
% legend('Diri','Diri Vel MG','Diri Pre MG',  'x-Peri', 'xy-Peri', 'Peri', 'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(10, -14.93, 'total no. of MG Sweeps');
% axis square 

% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('P2, 1-Vcycle/it','P3, 1-Vcycle/it', '-P2, 1-Vcycle/it', '-P3, 1-Vcycle/it', 'P2, 25-Vcycle/it','P3, 25-Vcycle/it', '-P2, 25-Vcycle/it', '-P3, 25-Vcycle/it','Location','NorthEast')
 text(25, -15.23, 'total no. GMRES It');
axis square 
%legend('Pure Dirichlet 1 MG/it', 'x- Periodic 1 MG/it', 'x,y -Periodic 1 MG/it', 'Pure Periodic 1 MG/it', 'Pure Periodic 2 MG/it', 'Pure Periodic 4 MG/it', 'Pure Periodic 10 MG/it', 'Location','NorthEastOutside')
%text(10, -12.43, 'total no. of MG Sweeps');

 hold on;
 text(-10.5, -7.15, 'log(r_n/r_0)');
 
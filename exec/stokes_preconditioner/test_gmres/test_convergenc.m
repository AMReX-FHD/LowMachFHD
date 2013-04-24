function test_convergenc

resid_1=load('restart_50.dat');
resid_2=load('restart_20.dat');
resid_3=load('restart_10.dat');

% periodic case 1-4 Jacobi iterataion, 5-8, GS iteration
plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b^'); 
hold on;

plot([0;resid_2(:,1)], [0;log10(resid_2(:,2))],  'r^');
hold on;

plot([0;resid_3(:,1)], [0;log10(resid_3(:,2))],  'k^');
hold on;


legend('Restart=50', 'Restart=20', 'Restart=10', 'Location','NorthEast');
text(220, -16.98, 'total no. of MG Sweeps');
clear all; format compact; format short

u_av=1.0; % Average concentration
dx=10.0; % Grid size (cell volume)

u_std = sqrt(u_av*dx)/dx; % Expected standard deviation

u_std_int = nearest(sqrt(u_av*dx))/dx; % Ensure that the bins are centered on actual integers
range = [u_av - 4*u_std_int, u_av + 6*u_std_int];
Nbins=nearest((range(2)-range(1))*dx)+1

bins = linspace(range(1), range(2), Nbins);
bins_th = linspace(range(1), range(2), 10000);
du = bins(2)-bins(1);
du_th = bins_th(2)-bins_th(1); % Compute the theory on a finer grid to get closer to zero
bounds=[bins(1), bins(end)]

Gaussian = exp(-(bins_th-u_av).^2/(2*u_std.^2))/sqrt(2*pi)/u_std;

posbins=bins_th(bins_th>0);
Poisson=zeros(1,length(bins_th));
Poisson(bins_th>0) = exp( -posbins.*(log(posbins/(u_av-1/(2*dx)))-1)*dx );
Poisson = Poisson / (sum(Poisson(:))*du_th);

compare=1 % Many curves at once or just a single run
if(compare) % Compare different
   
   f_data=load('/home/donev/HPC/FluctHydro/Chemistry/ReactDiff/Schlogl/fort.10.chem','-ascii'); f_data=f_data(:,2);
   u_hist_0=hist(f_data,bins); u_hist_0=u_hist_0/(sum(u_hist_0(:))*du);
   f_data=load('/home/donev/HPC/FluctHydro/Chemistry/ReactDiff/Schlogl/fort.10.diff','-ascii'); f_data=f_data(:,2);
   u_hist_1=hist(f_data,bins); u_hist_1=u_hist_1/(sum(u_hist_1(:))*du);
   f_data=load('/home/donev/HPC/FluctHydro/Chemistry/ReactDiff/Schlogl/fort.10.diff_chem','-ascii'); f_data=f_data(:,2);
   u_hist_2=hist(f_data,bins); u_hist_2=u_hist_2/(sum(u_hist_2(:))*du);
   
else

   f_data=load('/home/donev/HPC/FluctHydro/Chemistry/ReactDiff/Schlogl/fort.10','-ascii');
   f_data=f_data(:,2); % Remove the time label
   u_hist=hist(f_data,bins);
   u_hist_std=sqrt(u_hist);
   Z=1/(sum(u_hist(:))*du);
   u_hist = Z * u_hist ;
   u_hist_std = Z * u_hist_std ;
   
   u_mean_num = mean(f_data)
   u_std_num = std(f_data)
   ratio = u_std_num / u_std
   
end

figure(1); clf
if(0)
   errorbar(bins, u_hist, u_hist_std, 'ok'); hold on;
   plot(bins_th, Gaussian, '--r'); hold on;
   plot(bins_th, Poisson, '-g'); hold on;
   legend('Numerics','Gaussian','Poisson');
elseif(compare)
   semilogy(bins, u_hist_0, 'o--k'); hold on;
   semilogy(bins, u_hist_2, 's--b'); hold on;
   semilogy(bins, u_hist_1, 'd--m'); hold on;
   semilogy(bins_th, Gaussian, '-r'); hold on;
   semilogy(bins_th, Poisson, '-g'); hold on;
   %legend('Trap','Mid','CN','Gaussian','Poisson');
   legend('Chem only','Diff only','Chem+Diff','Gaussian','Poisson');  
else
   semilogy(bins, u_hist, 'ok'); hold on;
   semilogy(bins_th, Gaussian, '-r'); hold on;
   semilogy(bins_th, Poisson, '-g'); hold on;
   legend('Numerics','Gaussian','Poisson');
end
title('Monomodal Schlogl model N=10 particles per cell'); 

if(0) % Save the histograms to a file
array=zeros(Nbins,3);
array(:,1)=bins;
array(:,2)=Poisson;
array(:,3)=Gaussian;
save('/home/donev/HPC/FluctHydro/Chemistry/ReactDiff/hist-n_25.dat','array','-ascii');
end

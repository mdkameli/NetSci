close all
clear all
clc

% IMPORT & Polishing Data%
G = csvread('Adj.csv',1,1);
N = max(size(G));
A = sparse(G);
clear G;
figure(1)
spy(A)

%% pre-processing

Au = 1*(A+A'>0); % undirected network
Au = Au - diag(diag(Au)); % clear diagonal (you never know)

% Remove nodes which are NOT connected
pos = find(sum(Au)~=0);
A = A(pos,pos);
Au = Au(pos,pos);
spy(Au);

%% %%%%%%%%%%%%%%%%% EXTRACT THE DISTRIBUTION %%%%%%%%%%%%%%%%%%%%%%%%%

% Distribution
N = size(A,1);              % Number of Nodes
d = full(sum(A));           % Degree Vector
%d = d(d>0);                 % avoid zero degrees
Links_num = sum(d);         % Total number of links
k = unique(d);              % Degree Samples
pk = histc(d,k)';           % counts occurrences
pk = pk/sum(pk);            % normalize to 1

% Cumulative distribution
Pk = cumsum(pk,'reverse');

% Log Binning
klog = 10.^(0:0.1:ceil(log10(max(k))));
pklog = histc(d,klog)';     % counts occurrences
pklog = pklog/sum(pklog);   % normalize to 1

%% %%%%%%%%%%%%%%%%% MOMENTS OF DEGREE DISTRIBUTION %%%%%%%%%%%%%%%%%%%%

Mean_K = mean(k);           % First Moment of prob. distribution
Var_K = var(k);             % Second Moment of Prob. distribution (Express the spread)
Skew_K = skewness(k);       % Third Moment of Prob. distribution (How symmetric around average)


%% %%%%%%%%%%%%%%%%%%%%%%%% SHOW THE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure(2)
subplot(2,2,1)
plot(k,pk,'.')
grid
xlabel('k')
ylabel('PDF')
title('linear PDF plot')
subplot(2,2,2)
loglog(k,pk,'.')
grid
xlabel('k')
ylabel('PDF')
title('logarithmic PDF plot')
subplot(2,2,3)
loglog(klog,pklog,'.')
grid
xlabel('k')
ylabel('PDF')
title('logarithmic PDF plot (log bins)')
subplot(2,2,4)
loglog(k,Pk,'.')
grid
xlabel('k')
ylabel('CCDF')
title('logarithmic CCDF plot')


%% %%%%%%%%%%%%%%%%% PURE ML FITTING %%%%%%%%%%%%%%%%%%%%%%%%%

kmin = 40;
d2 = d(d>=kmin);                 % restrict range
ga = 1+1/mean(log(d2/kmin));     % estimate the exponent
disp(['gamma ML = ' num2str(ga)])


%% %%%%%%%%%%%%%%%%% SHOW THE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
subplot(2,2,2)
hold on
s1 = k.^(-ga); % build the PDF signal
loglog(k,s1/s1(103)*pk(103));
hold off
axis([xlim min(pk/2) 2*max(pk)])
subplot(2,2,3)
hold on
s1 = klog.^(-ga); % build the PDF signal
loglog(klog,s1/s1(21)*pklog(21));
hold off
axis([xlim min(pklog/2) 2*max(pklog)])
subplot(2,2,4)
hold on
s1 = k.^(1-ga); % build the CCDF signal
loglog(k,s1/s1(100)*Pk(100));
hold off
axis([xlim min(Pk/2) 2])


%% %%%%%%%%%%%%%%%%% ML FITTING WITH SATURATION %%%%%%%%%%%%%%%%%%
d1=d(d>30);
for ks = 1:max(k)
    kmin = min(d1);
    tmp = mean(log((d1+ks)/(kmin+ks)));
    ga2(ks) = 1+1/tmp;
    de(ks) = log(ga2(ks)-1)-log(kmin+ks)-ga2(ks)*tmp;
end
[~,ks] = max(de);
disp(['k_sat ML sat = ' num2str(ks)])
disp(['gamma ML sat = ' num2str(ga2(ks))])


%% %%%%%%%%%%%%%%%%% SHOW THE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
semilogy(de)
grid
xlabel('k_sat')
ylabel('ML target function')
title('best k_{sat} value')

figure(5)
% data
loglog(k,Pk,'.')
hold on
% ML fitting (we make sure that the plot follows the data)
s1 = k.^(1-ga); % build the CCDF signal
loglog(k,s1/s1(10)*Pk(10));
% ML fitting with saturation
s1 = ((k+ks)/(kmin+ks)).^(1-ga2(ks));
loglog(k,s1)
hold off
axis([xlim min(Pk/2) 2])
grid
xlabel('k')
ylabel('CCDF')
title('ML fittings')
legend('data','ML','ML with sat.')


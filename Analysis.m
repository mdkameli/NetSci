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
 display("Number of Nodes= "+N)
d = full(sum(A));           % Degree Vector
%d = d(d>0);                % avoid zero degrees
Links_num = sum(d);         % Total number of links --G (Let's decide if #node is this or we want to divide by 2)
 display("Total number of links= "+Links_num)

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

Mean_D = mean(d);           % First Moment of prob. distribution
                           % -- G before was Mean_K = mean(k), but actually k is
                           % -- G the vector of unique values, we should
                           % -- G use d; the same for all other mesures
display("First Moment of prob. distribution= "+Mean_D)


Var_D = var(d);             % Second Moment of Prob. distribution (Express the spread)

display("Second Moment of Prob. distribution= "+Var_D)

Skew_D = skewness(d);       % Third Moment of Prob. distribution (How symmetric around average)
display("Third Moment of Prob. distribution= "+Skew_D)


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
loglog(k,s1/s1(40)*Pk(40));
% ML fitting with saturation
s1 = ((k+ks)/(kmin+ks)).^(1-ga2(ks));
loglog(k,s1*exp(-1.2))
hold off
axis([xlim min(Pk/2) 2])
grid
xlabel('k')
ylabel('CCDF')
title('ML fittings')
legend('data','ML','ML with sat.')


%% %%%%%%%%%%%%%%%%%%%%%%%%%% ASSORTATIVITY %%%%%%%%%%%%%%%%%%%%%%%%%

k_tmp = (A*d')./d'; % the average degree of neighbours

% extract averages for each value of k
u = unique(d');
for k = 1:length(u)
    k_nn(k) = mean(k_tmp(d'==u(k))); 
end

% do the linear fitting
p = polyfit(log(u'),log(k_nn),1);
disp(['Assortativity factor ' num2str(p(1))]) %assortativity factor


%% %%%%%%%%%%%%%%%%% SHOW RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
loglog(d,k_tmp,'g.');
hold on
loglog(u,exp(p(2)+log(u)*p(1)),'r-');
loglog(u,k_nn,'k.');
hold off
grid
xlabel('k')
ylabel('k_{nn}')
title('Assortativity of the Collaboration Network')


%% %%%%%%%%%%%%%%% CLUSTRING COEFFICIENT 1 %%%%%%%%%%%%%%%%%%%%%%
%Clustring Coefficient: Measure density of links in the neighbourhood

cn = diag(full((Au*triu(Au)*Au)));          % Number of triangles for each node
Ei = zeros(size(d));
Ei = cn(d>1).';                             % Number of edges between nodes of neighbourhood
%E_Max = 0.5*d(d>1).*(d(d>1)-1);             % Maximum Number of Edges (Pairs)
C1 = 2*Ei./d(d>1)./(d(d>1)-1);                              % Clustring Coefficient
%C1 = C1(C1>0);
Cave1 = sum(C1)/N;

%% %%%%%%%%%%%%%% CLUSTRING COEFFICIENT 2 %%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Au)
  str(i).child = find(Au(i,:)>0);
end
N = length(Au);
deg = sum(Au(:,:));

% find neighboring links

for i=1:N
    adjn = Au(str(i).child,str(i).child);
    E(i) = sum(sum(adjn))/2;
    if deg(i)==1 | deg(i)==0
        C(i) = 0;
    else
        C(i) = 2*E(i)/deg(i)/(deg(i)-1);
    end

end
%C = C(C>0);
C = full(C);
Cave2 = sum(C)/N;
 display("Average Clustering Coefficient = "+Cave2)

%% CLUSTRING COEFFICIENT PROB. DISTRIBUTION
s = unique(C);                              % Unique Occurrences
%k = unique(d);
Cpk = histc(C,s)';                          % counts occurrences
Cpk = Cpk/sum(Cpk);                         % normalize to 1

% Cumulative distribution
CPk = cumsum(Cpk,'reverse');
figure(6);
plot(s,CPk,'+');
grid
hold on
hline = refline([0 Cave1]);
hline.Color = 'r';
hline.LineWidth = 1;
hline.DisplayName = 'Average Clustering Coefficient';
hold off
xlabel('k')
ylabel('PDF')
title('Clustring Coefficient Distribution')
%% %%%%%%%%%%%%%%%%% ROBUSTNESS %%%%%%%%%%%%%%%%%%%%
%Robustness: if you knock out x% of nodes/edges, how many % survive ?

% Inhomogeneity Ratio: 
%First we should check the availability of Giant Component that Molly-Reed
%Criterian holds

%Random Network
Au_update = Au;
Rand_inhom_Ratio = mean(d.^2) / mean(d);                    %A randomly wired network has a giant component

for i = 1:N-1
    j = ceil((N-i)*rand)+1;
    Au_update(:,j) = [];
    Au_update(j,:) = [];
    d_update = full(sum(Au_update));
    
    mom2_k = mean(d_update.^2);
    mom1_k = mean(d_update);
    
    Rand_inhom_Ratio_update = mom2_k ./mom1_k;
    Rand_inhom_Ratio = [Rand_inhom_Ratio Rand_inhom_Ratio_update];
end

figure(7);
loglog(Rand_inhom_Ratio,'.');
grid
hold on
hline = refline([0 2]);
hline.Color = 'r';
hline.LineWidth = 1;
hold off
xlabel('k')
ylabel('<K^2> / <K> ')
title('Robustness for Random failure')


% Scale free Network (gama = 2.5565)
% gama = ga;
% if gama > 3
%     ScFree_inhom_Ratio = kmin.*(gama-2)./(gama-3);
% elseif (2<gama) && (gama<3)
%         ScFree_inhom_Ratio = kmin.*(gama-2)./(3-gama).*N^((3-gama)./(gama-1));
% elseif (1<gama) && (gama<2)
%         ScFree_inhom_Ratio = kmin.*(2-gama)./(3-gama).*N^(1./(gama-1));
% end
%display('For the Scale free Network The Inhomogeneity Ratio is = '+ScFree_inhom_Ratio);
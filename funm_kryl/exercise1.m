%% Initialization
clear all; 
clc;
close all;
%% Load of the graph
load Harvard500;
%% Visualize the sparsity structure of the graph
spy(Problem.A);
%% Create a directed graph
NodeName=cellstr(Problem.aux.nodename);
G=digraph(Problem.A',NodeName);

%% Compute Pagerank
pr = centrality(G,'pagerank','FollowProbability',0.85);

%% Assign to nodes the computed PageRank 
G.Nodes.PageRank = pr;

%% Assign to nodes Indegree and OutDegree

G.Nodes.InDegree = indegree(G);
G.Nodes.OutDegree = outdegree(G);

%% Visualize Obtained Results

G.Nodes(1:25,:)
%% Extract and plot a subgraph containing all nodes whose score is greater than 0.01.
selection=find(G.Nodes.PageRank > 0.01);
% H = subgraph(G,nodeIDs) returns a subgraph of G that contains only the nodes specified by nodeIDs.
H = subgraph(G,selection);
figure(2)
% plot(H,'NodeLabel',nodename(selection),'NodeCData',H.Nodes.PageRank,'Layout','force');
plot(H,'NodeLabel',{},'NodeCData',H.Nodes.PageRank,'Layout','force');
colorbar

%% Visualize the results according PageRank
[spr,per]=sort(pr,'descend')
result= table;
result.NodeName=NodeName(per);
result.PageRank=spr;
result.InDegree = G.Nodes.InDegree(per);
result.OutDegree = G.Nodes.OutDegree(per);
result(1:25,:)


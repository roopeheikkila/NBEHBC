%In this pre-assignment the task is to test the Brain Connectivity Toolbox
%The macaque network from https://sites.google.com/site/bctnet/datasets was
clear;

%load the dataset
load('macaque71.mat');

%% Degree and Similarity

% Degree:
% degrees_und.m (BU, WU networks); degrees_dir.m (BD, WD networks).
degrees = degrees_und(CIJ);
figure(1);
subplot(3,1,1)
stem(degrees)
subplot(3,1,2)

histogram(degrees)
subplot(3,1,3)
boxplot(degrees)
% 
% Strength: 
% strengths_und.m (WU networks); strengths_dir.m (WD networks).
% strengths_und_sign.m (WU signed networks).

% Joint degree
% jdegree.m (BD, WD networks).


% Topological overlap: 
% gtom.m (BU networks).

% 
% Neighborhood overlap: 
% edge_nei_overlap_bu.m (BU networks); edge_nei_overlap_bd.m (BD networks).


% Matching index: 
% matching_ind_und.m (BU networks); matching_ind.m (BD networks).


%% Centrality

%Betweenness centrality:
bc =betweenness_bin(CIJ);
plot_measure(bc, 'betweenness centrality','betweenness centrality');

%Edge betweenness centrality
ebc = edge_betweenness_bin(CIJ);
figure();
%create a heatmap where each point is a edge between two nodes
imagesc(ebc);
ylabel('node index');
xlabel('node index');
title('edge betweenness centrality');

%Within-module degree z-score:
%TODO
%Participation and related coefficients:
%TODO

%Eigenvector centrality: 
ec = eigenvector_centrality_und(CIJ);
plot_measure(ec, 'eigenvector centrality', 'Eigenvector centrality');
%%
%PageRank centrality:
prc = pagerank_centrality(CIJ,0.5);
plot_measure(prc, 'pagerank centrality', 'Pagerank centrality');

%%
[x,y] = adjacency_plot_und(CIJ);
plot(x,y)

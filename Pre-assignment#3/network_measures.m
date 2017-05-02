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
bc =betweenness_bin(CIJ);figure(1);

subplot(3,1,1)
    stem(bc)
    xlabel('node index')
    ylabel('betweenness centrality');
    title('betweenness centrality');
subplot(3,1,2)
    histogram(bc);
    xlabel('betweenness centrality');
    ylabel('amount of nodes');
subplot(3,1,3)
    boxplot(bc)
    ylabel('betweenness centrality');


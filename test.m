% Author: David Martinez Rego
% Example of how to use the LSH routine for euclidean
% distances

clear all
close all

rng('shuffle')

addpath('./lib')
mu = randi(5, [1,2]);
sigma = diag(1*ones(1, 2));
rng('shuffle')

r = mvnrnd(mu,sigma,2000);
q = mvnrnd(mu,sigma,2);
plot(r(:, 1), r(:, 2), 'b+');
hold on
plot(q(:, 1), q(:, 2), 'r.');
pause
    
% Usage
% R : radious
% succP : probability of success
% queries : query points
% set : dataset
% maxReported: maximum number of reported NNs
%
% return : res --> list of (point, number of NNs for that point, NNs list), ...

% [res] = lshfind( R, succP, queries, set, maxReported);
a = lshfind(0.5, 0.9, q', r', 4)

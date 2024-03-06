 close all;
clear all;
clc
clear
%% parameter setting
format long
addpath('utils/');
addpath('metric_utils\');

%% begin to process one image sequence

readPath = '.\data';
savePath = '.\result';
p = 1;
frame = 3;
lambdaL = 30;
mu= 0.010;   
per = 1;

if ~exist('.\result')
    mkdir('.\result')
end

tuneopts.temporal_step = frame;
tuneopts.per = per;
tuneopts.lambdaL = lambdaL;
tuneopts.mu = mu;
target_detection(char(readPath), savePath, tuneopts);


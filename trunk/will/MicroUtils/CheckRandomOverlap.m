%function results=CheckRandomOverlap(NUM_CHECKS,varargin)
%
%
%
%
%
%
%
%
%
%
% 
%

NUM_CHECKS=2;
params=1:NUM_CHECKS;
result=startmulticoremaster(@CheckRandom,params,'c:\distcomp');
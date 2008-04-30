function [message,results,timeElapsed] = runcontest(varargin)
%RUNCONTEST Test a contest entry.
%   [MESSAGE,RESULTS,TIME] = RUNCONTEST runs the M-file solver.m against all
%   the problems defined in testsuite_sample.mat. MESSAGE returns a summary
%   of the testing.  RESULTS measures how well the entry solved the problem
%   and TIME measures the time the entry took to compute its answer.
%
%   ... = RUNCONTEST(TRUE) steps into each of the responses.

% The MATLAB Contest Team
% November 2007

load testsuite_sample testsuite
n = numel(testsuite);
responses = cell(n,1);
scores = zeros(n,1);

if nargin == 0 % no board drawing
    time0 = cputime;
    for k = 1:n
        inputs = struct2cell(testsuite(k));
        responses{k} = solver(inputs{:});
    end
    timeElapsed = cputime-time0;
    for k = 1:numel(testsuite)
        inputs = struct2cell(testsuite(k));
        scores(k) = grade(inputs{:},responses{k});
    end
else          % step into each of the responses
    for k = 1:n
        inputs = struct2cell(testsuite(k));
        responses{k} = solver(inputs{:});
        scores(k) = visualize(inputs{:},responses{k});
        pause
    end
    timeElapsed = NaN;
end

results = sum(scores);
message = sprintf('results: %.4f\ntime: %.2f',results,timeElapsed);

function [data,time_vect] = generate_data(A,seconds,sample_rate,state_timecourse)
%%function [data,time_vect] = generate_data(A,seconds,sample_rate,state_timecourse)
%
% Function generating data from an autogressive parameter matrix with
% optional switching between states.
%
% Inputs:
%   A: matrix
%       Autoregressive parameter matrix [nnodes x nnodes x order x state]
%       Assumes lag 1 is eye
%  seconds: float
%       Time window of data to generate
%  sample_rate: float
%       sampling frequency of the system
%  state_timecourse: vector
%       Optional vector allocating each time point to one of the states in
%       the final dimension of A. If ommitted we allocate everything to
%       state 1.
%
% Returns:
%   data: matrix
%       [nnodes, seconds*sample_rate] matrix of data
%   time_vect: vector
%       Time stamps for each sample in data.
%
% AQ

% Allocate everything to state 1 unless told otherwise
if nargin < 4 || isempty(state_timecourse)
    state_timecourse = ones(seconds*sample_rate,1);
end

% Make time-vector
time_vect = linspace(0,1/sample_rate,seconds*sample_rate);
nsamples = length(time_vect);

% Get some info from A
[~,nnodes,order,nstates] = size(A);

% Generate white noise
data = randn(nnodes,nsamples);

% Make the data
for idx = order:length(time_vect)
    for ilag = 2:order
        data(:,idx) = data(:,idx) - A(:,:,ilag,state_timecourse(idx))*data(:,idx-(ilag-1));
    end
end

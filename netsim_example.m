%% An example simulation
%
% Here we generate a three node network which transitions between three
% spectrally-defined states. Each mode is an oscillation with a distribution
% across the network. Each state has a different set of oscillatory modes with
% different amplitudes and relative phases aacross the three nodes in the
% network.
%
% A mode is defined as a peak frequency, global amplitude, amplitude-per-node
% and relative-phase-per-mode. These are defined below.
%
% mode_info defines a network within one mode.
%   mode_info{1}.freq: defines the peak frequency for the mode. A value of 0 indicates 1/f noise
%   mode_info{1}.mode_amp: defines the pole amplitude attached to the freq
%   mode_info{1}.node_amp: defines the relative amplitude of the oscillation for each node
%   mode_info{1}.phase: defines the relative phase lag of the oscillation for each node.
%
% The number of channels is implicitly defined by the number of values in
% node_amp and phase
%
% It is best to have the same number of mode and channels per state.
% Signal-to-noise can be manipulated by changing the node_amp values directly
%clear all; close all

% Define state 1
clearvars mode_info
mode_info{1}.freq = 0;
mode_info{1}.mode_amp = .9;
mode_info{1}.node_amp = [1 1 1];
mode_info{1}.phase = [ pi/3 0 0 ];
mode_info{2}.freq  = 10;
mode_info{2}.mode_amp   = .94;
mode_info{2}.node_amp   = [ .5 1.5 1 ];
mode_info{2}.phase = [ 0 pi/2 0 ];
mode_info{3}.freq  = 27;
mode_info{3}.mode_amp   = .9;
mode_info{3}.node_amp   = [ 1 1 1 ];
mode_info{3}.phase = [ 0 0 .1 ];
signal{1} = mode_info;

% Define state 2
clearvars mode_info
mode_info{1}.freq = 0;
mode_info{1}.mode_amp = .9;
mode_info{1}.node_amp = [1 1 1];
mode_info{1}.phase = [ pi 0 0 ];
mode_info{2}.freq  = 13;
mode_info{2}.mode_amp   = .94;
mode_info{2}.node_amp   = [ .2 .3 1 ];
mode_info{2}.phase = [ 0 0 pi ];
mode_info{3}.freq  = 33;
mode_info{3}.mode_amp   = .9;
mode_info{3}.node_amp   = [ 1 1 1 ];
mode_info{3}.phase = [ 0 0 0 ];
signal{2} = mode_info;

% Define state 3
clearvars mode_info
mode_info{1}.freq = 0;
mode_info{1}.mode_amp = .9;
mode_info{1}.node_amp = [1 1 1];
mode_info{1}.phase = [ pi 0 0 ];
mode_info{2}.freq  = 10;
mode_info{2}.mode_amp   = .94;
mode_info{2}.node_amp   = [ .2 1.5 1 ];
mode_info{2}.phase = [ 0 pi 0 ];
mode_info{3}.freq  = 40;
mode_info{3}.mode_amp   = .9;
mode_info{3}.node_amp   = [ 0 1 1 ];
mode_info{3}.phase = [ 0 0 0 ];
signal{3} = mode_info;

sample_rate = 128;
% Generate modal parameters and short data for estimating spectra
[A(:,:,:,1),~,data1] = generate_modal_network(signal{1},sample_rate,100);
[A(:,:,:,2),~,data2] = generate_modal_network(signal{2},sample_rate,100);
[A(:,:,:,3),~,data3] = generate_modal_network(signal{3},sample_rate,100);

% Plot state-wise power spectra
figure;
subplot(311);hold on; grid on
[pxx,f] = pwelch(data1(1,:),[],[],[],sample_rate);plot(f,20*log(pxx));
[pxx,f] = pwelch(data1(2,:),[],[],[],sample_rate);plot(f,20*log(pxx));
[pxx,f] = pwelch(data1(3,:),[],[],[],sample_rate);plot(f,20*log(pxx));
title('State 1');

subplot(312);hold on; grid on
[pxx,f] = pwelch(data2(1,:),[],[],[],sample_rate);plot(f,20*log(pxx));
[pxx,f] = pwelch(data2(2,:),[],[],[],sample_rate);plot(f,20*log(pxx));
[pxx,f] = pwelch(data2(3,:),[],[],[],sample_rate);plot(f,20*log(pxx));
title('State 2');

subplot(313);hold on; grid on
[pxx,f] = pwelch(data3(1,:),[],[],[],sample_rate);plot(f,20*log(pxx));
[pxx,f] = pwelch(data3(2,:),[],[],[],sample_rate);plot(f,20*log(pxx));
[pxx,f] = pwelch(data3(3,:),[],[],[],sample_rate);plot(f,20*log(pxx));
title('State 3');
xlabel('Frequency (Hz)');
ylabel('dB');

legend({'Node 1','Node 2','Node 3'});



%% Generate mixed data

% Create the state lifetimes, drawn from a gamma distribution.
state_per{1} = fix(gamrnd(5,8,10000,1));  % Swich approx every 250ms
state_per{2} = fix(gamrnd(30,8,1000,1)); % Switch approx every second
state_per{3} = fix(gamrnd(100,8,1000,1)); % Switch approx every three seconds

% The order of the state occurances
state_ts_val = randi([1 3],10000,1 ); % state allocations

% Data parameters
seconds = 600; % make 10 minutes of data

% Exapand lifetimes and state allocs into full state-timecourse
time_vect = linspace(0,1/sample_rate,seconds*sample_rate);
state_ts{1} = ones(seconds*sample_rate,1);
state_ts{2} = ones(seconds*sample_rate,1);
state_ts{3} = ones(seconds*sample_rate,1);

for idx = 1:3

    start = 1;
    state = 1;

    while start+state_per{idx}(state) < seconds*sample_rate

        state_ts{idx}(start:start+state_per{idx}(state)) = state_ts_val(state);

        start = start + state_per{idx}(state);
        state = state + 1;

    end

end

% Generate the mixed data at each time-scale
clearvars data
data{1} = generate_data(A,seconds,sample_rate,state_ts{1});
data{2} = generate_data(A,seconds,sample_rate,state_ts{2});
data{3} = generate_data(A,seconds,sample_rate,state_ts{3});

save('Modal_simulation','data','sample_rate','time_vect','signal','state_ts');

%% Plot

% Expand state time-course into easily plottable form
for iscale = 1:3
    statewise{iscale} = zeros(3,seconds*sample_rate);
    for istate = 1:3ls
        
        statewise{iscale}(istate,:) = state_ts{iscale}' == istate;
    end
end

% Don't plot everything...
stop = 50*sample_rate;

order = size(A,3);

figure;
subplot(611);hold on; grid on;
plot(time_vect(1:stop),data{1}(:,1:stop)');axis('tight');
title('Fast Time-Series');
subplot(612);
area(time_vect(order:stop+order),statewise{1}(:,order:stop+order)');axis('tight');
title('Fast States');

subplot(613);hold on; grid on;
plot(time_vect(1:stop),data{2}(:,1:stop)');axis('tight');
title('Medium Time-Series');
subplot(614);
area(time_vect(order:stop+order),statewise{2}(:,order:stop+order)');axis('tight');
title('Medium States');

subplot(615);hold on; grid on;
plot(time_vect(1:stop),data{3}(:,1:stop)');axis('tight');
title('Slow Time-Series');
subplot(616);
area(time_vect(order:stop+order),statewise{3}(:,order:stop+order)');axis('tight');
title('Slow States');


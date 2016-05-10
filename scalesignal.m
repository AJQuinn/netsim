function [signal,rescale_db,accuracy] = scalesignal(signal, diff_db, ref_signal);
%function [signal,rescale_db,accuracy] = scalesignal(signal, diff_db, ref_signal);
%
% Scale the variance of a signal to a set difference decibels relative to a
% reference signal

    % Get original signal power
    sig_db = 20 * log10( sqrt ( sum ( power(signal,2) ) ) );

    % Calculate actual dB scaling if passed a reference signal
    if nargin == 3
        ref_db = 20 * log10( sqrt ( sum ( power(ref_signal,2) ) ) );
        diff_db = ref_db - sig_db + diff_db;
    end

    % Calculate scaling factor
    scale_factor = 10 ^ ( diff_db / 20);

    % Scale signal
    signal = scale_factor * signal;

    % Calculate scaled signal db
    rescale_db = 20 * log10( sqrt ( sum ( power(signal,2) ) ) );

    % Estimate accuracy
    accuracy = (sig_db - diff_db) / rescale_db;

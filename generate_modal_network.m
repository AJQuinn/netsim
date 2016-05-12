function [ A, modal_matrix, data ] = generate_modal_network( mode_info, sample_rate, seconds, noise_db, verbose )
%%function [ A, modal_matrix, data ] = generate_modal_network( mode_info, seconds, sample_rate )
%
% Generates the autoregressive parameters for a system based on a set of
% input modes and a mixing matrix. Can optionally generate data if
% requested.
%
% Inputs:
%   mode_info: Cell array of structs
%        Cell array containing the modal definitions
%        mode_info{n}: info for nth mode
%        mode_info{n}.freq: list of frequencies (ie [5 18]) bounded between
%                            0 and nyquist
%        mode_info{n}.mode_amp: scalar value < 1 indicating the mode pole magnitude
%        mode_info{n}.node_amp: list of amplitudes of corresponding freqs (ie [.9
%                         .5]) bounded between 0 and 1
%        mode_info{n}.phase: list of relative phase differences in rads for
%                            corresponding freqs (ie [0 2.5]) bounded between
%                            -pi and pi. Default is 0.
%   sample_rate: float
%        Sampling frequency of the system in Hz
%   seconds: float
%        How much data to generate, if 0 or empty, nothing is generated.
%   noise_db: int
%        Noise level in decibels relative to signal, if empty, no noise is
%        added
%
% Returns:
%   A: matrix
%       Autoregressive parameter matrix [nnodes x nnodes x modelorder]
%   modal_matrix: matrix
%       matrix containing the relative phases within the system
%   data: matrix
%       Generated data, is seconds > 0
%
% AQ

% Housekeeping
nnodes = length(mode_info{1}.node_amp);
nmodes = length(mode_info);
nyq = sample_rate/2;
mode_poles = [];
if nargin < 5; verbose = false; end

% Generate the poles and phases
for imode = 1:length(mode_info)

    % Frequency as proportion of nyqyuist
    mode_theta = (mode_info{imode}.freq/nyq) * pi;
    % Pole location in z-space
    mode_poles = [mode_poles ...
                  mode_info{imode}.mode_amp * (cos(mode_theta) + 1i*sin(mode_theta)) ];
    % Pole complex conjugate if the pole isn't strictly real
    if mode_info{imode}.freq > 0
        mode_poles = [mode_poles ...
                  mode_info{imode}.mode_amp * (cos(mode_theta) - 1i*sin(mode_theta)) ];
    end

    % Complex value indicating relative amplitude and phase per node
    modal_v = mode_info{imode}.node_amp' .* ...
                (cos(mode_info{imode}.phase) + 1j*sin(mode_info{imode}.phase))';

    % Stack full modal matrix
    if imode == 1
        modal_matrix = modal_v;
    else
        modal_matrix = cat(2,modal_matrix,modal_v);
    end

    if mode_info{imode}.freq > 0
        modal_matrix = cat(2,modal_matrix,conj(modal_v));
    end

end

C_valid = false;
while C_valid == false

    % Add some noise if there are fewer than mp poles
    % this is to ensure that the eigenvectors will have a well formed inverse
    % This is noisy, so we may repeat if the time-parameters are not
    % acceptable.
    npoles = size(modal_matrix,2);
    if npoles < nmodes*nnodes
        modal_matrix_full = cat(2,modal_matrix,randn(nnodes,(nnodes*nmodes)-length(mode_poles))/25);
        mode_poles_full = cat(2,mode_poles,randn(1,(nnodes*nmodes)-length(mode_poles))/25);
    end

    % Create eigenvector matrix in Vandermonde form
    eigenvecs = zeros(nmodes*nnodes,nmodes*nnodes);
    eigenvecs((nmodes*nnodes)-nnodes+1:nmodes*nnodes,:) = modal_matrix_full;
    if length(mode_info) > 1
        for idx = 1:length(mode_info)-1
            for ipole = 1:nnodes*nmodes
                eigenvecs(((nmodes-idx)*nnodes)-nnodes+1:(nmodes-idx)*nnodes,ipole) = ...
                    modal_matrix_full(:,ipole).*(mode_poles_full(ipole).^idx);
            end
        end
    end

    % Create time domain parameters in companion form
    C = eigenvecs*diag(mode_poles_full)*pinv(eigenvecs);

    % Sanity check, C should be real to reasonable precision
    if sum(sum(imag(C))) > 1e-13;
        disp(sum(sum(imag(C))))
        test1 = false;
        if verbose==true;warning('Parameters have large imaginary component!');end
    else
        test1 = true;
    end

    C = real(C);

    % Sanity check, The bottom rows of C should be sparse ones
    if abs( sum(sum(C(nnodes+1:end,:))) - ( nnodes * (nmodes - 1) ) ) >  1e-13
        test2 = false;
        if verbose==true;warning('Parameters are not in full companion form!');end
    else
        test2 = true;
    end

    % Accept parameters if we pass both sanity checks
    if test1 == true && test2 == true
        C_valid = true;
    end

end

% Companion to standard form
A = A_form_swap(-C,nnodes,size(C,1)/nnodes);

% Generate dataset if requested
if seconds > 0
    data = generate_data(A,seconds,sample_rate);

    if nargin == 4 && ~isempty(noise_db)
        for idx = 1:nnodes
            data(idx,:) = data(idx,:) + scalesignal( randn(1,size(data,2)),noise_db,data(idx,:) );
        end
    end

end

end


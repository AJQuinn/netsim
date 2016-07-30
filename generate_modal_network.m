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
nmodes = 0;
nyq = sample_rate/2;
mode_poles = [];
if nargin < 5; verbose = true; end

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
    if isreal(mode_info{imode}.node_amp)
        modal_v = mode_info{imode}.node_amp' .* ...
                (cos(mode_info{imode}.phase) + 1j*sin(mode_info{imode}.phase))';
    else
        modal_v = mode_info{imode}.node_amp.';
    end

    % Stack full modal matrix
    if imode == 1
        modal_matrix = modal_v;
    else
        modal_matrix = cat(2,modal_matrix,modal_v);
    end

    if mode_info{imode}.freq > 0
        modal_matrix = cat(2,modal_matrix,conj(modal_v));
        nmodes = nmodes + 2;
    else
        modal_matrix(:,end) = abs(modal_matrix(:,end));
        nmodes = nmodes + 1;
    end

end

% Add some noise if there are fewer than mp poles
if nnodes == 1
    % Single channel, order equals modes
    p = nmodes;
elseif nmodes < nnodes
    % Fewer modes than nodes
    p = (floor(nmodes/nnodes)) + 1;
else
    % More modes than nodes
    p = (floor(nmodes/nnodes));
end

poles_to_add = (nnodes*p)-nmodes;
if poles_to_add > 0

    disp('Adding dummy poles, may have unexpected results');

    % Create filler poles, evenly spaced across the rest of the
    % spectrum.
    start_angle = max(angle(mode_poles));
    % Find angles for new filler poles, relative to largest frequency
    % defined poles
    new_angles = linspace(start_angle, 2*pi - start_angle, poles_to_add+2 );
    new_angles = new_angles(2:end-1); % remove defined poles
    % create filler pole coordinates at fixed magnitude
    new_poles = .25*( cos(new_angles) + 1j*sin(new_angles) );
    % Sort into increasing frequency, should make conjugate pair adjacent
    [~,I] = sort(abs(angle(new_poles)));
    new_poles = new_poles(I);
    % Add them to the list
    mode_poles_full = cat(2,mode_poles,new_poles);

    % Generate orthogonal vectors to use as eigenvectors
    imbasis = orth(reshape( (1:prod(poles_to_add*nnodes))-nnodes, nnodes, poles_to_add) );
    rebasis = orth(reshape( (prod(poles_to_add*nnodes):-1:1)-nnodes, nnodes, poles_to_add) );
    imbais = gallery('orthog',nnodes,-2);
    rebasis = gallery('orthog',nnodes,-2)';

    % Create filler eigenvectors, very small weighting into the spatial
    % dimensions. Take care of new conjugate pairs here
    modal_matrix_full = modal_matrix;
    for idx = 1:floor(poles_to_add/2)
        % add matching vecs for conjugate poles
        new_vec = complex(rebasis(:,idx),imbasis(:,idx)) / 1000;
        modal_matrix_full = cat(2,modal_matrix_full,...
                [new_vec conj(new_vec)]);
    end
    % Add vectors for a real pole if necessary
    if mod(poles_to_add,2) == 1
    if imag(new_poles(end)) < 1e-10
        modal_matrix_full = cat(2,modal_matrix_full,complex(rebasis(:,end),imbasis(:,end))/1000);
    end
    end

else

    modal_matrix_full = modal_matrix;
    mode_poles_full = mode_poles;

end
nmodes = length(modal_matrix_full);

% Create eigenvector matrix in Vandermonde form
eigenvecs = zeros(nmodes,nmodes);
eigenvecs((nmodes)-nnodes+1:nmodes,:) = modal_matrix_full;
if nmodes > nnodes
    %for idx = 1:( (nmodes / nnodes) - 1 )
    for idx = ( (nmodes / nnodes) - 1 ):-1:1
        for ipole = 1:length(mode_poles_full)
            eigenvecs((1:nnodes)+((idx-1)*nnodes),ipole) = ...
                modal_matrix_full(:,ipole).*(mode_poles_full(ipole).^idx);
        end
    end
end

% Sanity check, eigenvectors should be linearly independent
if rank(eigenvecs) < size(eigenvecs,1)
    [~,jb] = rref(eigenvecs);
    disp(['Eigenvectors are low rank:' num2str(rank(eigenvecs)) ' < ' num2str(size(eigenvecs,1))])
    disp(['Rows ' num2str(setdiff(1:nmodes,jb)) ' potentially co-linear'])
    error('Eigenvectors are low rank')
end

% Sanity check, is the eigenvector matrix properly invertible
% we could also just check cond, but the absolute value at which the
% condition is 'high' is hard to say. Better to absolutely determine the
% error
if abs(sum(sum(eigenvecs*inv(eigenvecs) - eye(nmodes)))) > 1e-10
    disp(['Right Eigenvector matrix is poorly conditioned: ' num2str(cond(eigenvecs))]);
end

left = inv(eigenvecs');
if abs(sum(sum(left*inv(left) - eye(nmodes)))) > 1e-10
    disp(['Left Eigenvector matrix is poorly conditioned: ' num2str(cond(eigenvecs))]);
end

% Create time domain parameters in companion form
C = eigenvecs*diag(mode_poles_full)*inv(eigenvecs);
%mode_poles_full

% Sanity check, can we reconstruct with left eigenvectors
Cl = inv(left') * diag(mode_poles_full) * left';
if sum(sum(abs(C - Cl))) > 1e-10
    disp('Left eigenvectors do not recover companion form!');
end

% Sanity check, C should be real to reasonable precision
if sum(sum(abs(imag(C)))) > 1e-10;
    disp(sum(sum(imag(C))))
    test1 = false;
    if verbose==true;warning('Parameters have large imaginary component!');end
else
    test1 = true;
end

C = real(C);

% Sanity check, The bottom rows of C should be sparse ones
if  sum(sum(abs(C(nnodes+1:end,:)))) - ( nmodes - nnodes ) > 1e-10;
    test2 = false;
    if verbose==true;warning('Parameters are not in full companion form!');end
else
    test2 = true;
end

% Accept parameters if we pass both sanity checks
if test1 == true && test2 == true
    C_valid = true;
else
    error('Generated parameter matrix is not valid, please check input options');
end

% Companion to standard form
A = A_form_swap(-C,'comp2full',nnodes,size(C,1)/nnodes);

% Generate dataset if requested
if seconds > 0
    data = generate_data(A,seconds,sample_rate);

    if nargin >= 4 && ~isempty(noise_db)
        for idx = 1:nnodes
            data(idx,:) = data(idx,:) + scalesignal( randn(1,size(data,2)),noise_db,data(idx,:) );
        end
    end

end

end


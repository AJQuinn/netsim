function [ Anew ] = A_form_swap( A, nnodes, order )
%% function [ Anew ] = A_form_swap( A, nnodes, order )
%
% Convert autoregressive parameter matrix between 3D typical form and 2D
% companion form
%
% Input
%   A - autoregressive parameters
%   nnodes - number of time-series
%   order - number of time lags in model


if ndims(A) == 2

    % Change to 3d form
    Anew = reshape(A(1:nnodes,:),nnodes,nnodes,[]);
    Anew = cat(3,eye(nnodes),Anew);

else

    % Change to companion
    Anew = zeros(nnodes*(order),nnodes*(order));
    Anew(1:nnodes,:) = reshape(A(:,:,2:end),nnodes,[]);
    Anew(nnodes+1:end,1:end-nnodes) = eye(-nnodes + (order)*nnodes);

end

end


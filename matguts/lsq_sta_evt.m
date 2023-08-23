function [ sta_terms,evt_terms ] = lsq_sta_evt( all_mat, epsilonS,epsilonE )
% [ sta_terms,evt_terms ] = lsq_sta_evt( all_mat, epsilonS,epsilonE )
% Script to solve the damped least squares problem for best fitting station and
% event values from some matrix of values that are (at each station) a
% combination of station and event terms.
% 
%  Usage:
%       [ sta_terms,evt_terms ] = lsq_sta_evt( all_mat, epsilonS,epsilonE  )
% 
%  Inputs:
%       all_mat   = nstas x nevts matrix of values (possibly with nans)
%       epsilonS   = strength of damping for station terms.
%       epsilonE   = strength of damping for event terms.
% 
%  Outputs:
%       sta_terms = nstas x 1 vector of best-fitting station values
%       evt_terms = nevts x 1 vector of best-fitting event values

if nargin < 2
    epsilonS = 0.00001;
end
if nargin < 3
    epsilonE = 0.00001;
end

[nstas,nevts] = size(all_mat);

N = sum(sum(~isnan(all_mat))); % data
M = nstas+nevts;               % model parameters


si = zeros(2*N,1);
sj = zeros(2*N,1);
s = zeros(2*N,1);

d = zeros(N,1);

kk = 0;
for ie = 1:nevts
    for is = 1:nstas
        if isnan(all_mat(is,ie)), continue, end
        kk = kk+1; % n_obs
        si(2*kk+[-1 0]) = kk;
        sj(2*kk+[-1 0]) = [is, nstas+ie];
        s(2*kk+[-1 0]) = 1;
        
        d(kk) = all_mat(is,ie);
    end
end

% 'data' kernel relating each val to the sum of a station and event term
G = sparse(si,sj,s,N,M,2*N);

% damping matrix
L = [epsilonS*sparse(1:nstas,1:nstas,ones(nstas,1),nstas,nstas),sparse(nstas,nevts);
     sparse(nevts,nstas),epsilonE*sparse(1:nevts,1:nevts,ones(nevts,1),nevts,nevts)];

% constraint line - station terms must average to zero;
C = sparse([ones(1,nstas),zeros(1,nevts)]);

F = [G; L; C];
f = [d; zeros(nevts+nstas,1);0];

m = [F'*F]\F'*f;

sta_terms = m(1:nstas);
evt_terms = m([1:nevts] + nstas);

% E = (d - G*m)'*(d - G*m)/N

% [sta_terms,nanmean(all_mat,2)]

end
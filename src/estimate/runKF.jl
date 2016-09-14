#  Replication files for:
#  ""Nowcasting", 2010, (by Marta Banbura, Domenico Giannone and Lucrezia Reichlin),
#  in Michael P. Clements and David F. Hendry, editors, Oxford Handbook on Economic Forecasting.
#
#  The software can be freely used in applications.
#  Users are kindly requested to add acknowledgements to published work and
#  to cite the above reference in any resulting publications
# #--------------------------------------------------------------------------
# # KALMAN FILTER
# #--------------------------------------------------------------------------
"""
```
KalmanType{T} <: AbstractModel{T}
```
        S.Am       Predicted state vector  A_t|t-1      (nobs x m)
        S.AmU      Filtered  state vector  A_t|t        (nobs+1 x m)
        S.Pm       Predicted covariance of A_t|t-1      (nobs x m x m)
        S.PmU      Filtered  covariance of A_t|t        (nobs+1 x m x m)
        S.loglik   Value of likelihood function

        S.Am       Estimates     a_t|t-1                  (nobs x m)
        S.Pm       P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
        S.AmU      Estimates     a_t|t                    (nobs x m)
        S.PmU      P_t|t   = Cov(a_t|t)               (nobs x m x m)
        S.AmT      Estimates     a_t|T                    (nobs x m)
        S.PmT      P_t|T   = Cov(a_t|T)               (nobs x m x m)
        S.PmT_1    Cov(a_ta_t-1|T)
        where m is the dim of state vector and t = 1 ...T is time
"""

type KalmanType{T} # NOT of the abstract type Model... X !<: AbstractModel{T}
      Am::Array{Float64,2}       # Predicted state vector  A_t|t-1 (nobs   x m)
      AmU::Array{Float64,2}      # Filtered  state vector  A_t|t   (nobs+1 x m)
      Pm::Array{Float64,2}       # Predicted covariance of A_t|t-1 (nobs   x m x m)
      PmU::Array{Float64,2}      # Filtered  covariance of A_t|t   (nobs+1 x m x m)
      AmT::Array{Float64,2}      # Estimates     a_t|T             (nobs   x m)
      PmT::Array{Float64,2}      # P_t|T   = Cov(a_t|T)            (nobs   x m x m)
      PmT_1::Array{Float64,2}    # Cov(a_t'*a_t-1|T)               (nobs   x m x m)
      loglik::Float64            # Value of likelihood function
  end

#TODO need to define all the input types to the functions
# will need to initialize a S of KalmanType with empty values of the correct size if perfect
# REWRITE THE function definition to take a KalmanType and mutate it
function runKF!(S::KalmanType, y::Array{Float64,2},
                A::Array{Float64,2}, C::Array{Float64,2},
                Q::Array{Float64,2}, R::Array{Float64,2},
                x_0::Array{Float64,2}, Sig_0::Array{Float64,2});

SKF!(S,y,C,R,A,Q, x_0, Sig_0); # mutates and outputs the S type
FIS!(S,y,C,R,A,Q);  # mutates and outputs the S type

xsmooth = S.AmT;
Vsmooth = S.PmT;
VVsmooth = S.PmT_1;
loglik = S.loglik;

return xsmooth, Vsmooth, VVsmooth, loglik
end

#______________________________________________________________________
function SKF!(S,Y,Z,R,T,Q,A_0,P_0)
#______________________________________________________________________
# Kalman filter for stationary systems with time-varying system matrices
# and missing data.
#
# The model is        y_t   = Z * a_t + eps_t
#                     a_t+1 = T * a_t + u_t
#
#______________________________________________________________________
# INPUT
#        Y         Data                                 (nobs x n)
# OUTPUT
#        S.Am       Predicted state vector  A_t|t-1      (nobs x m)
#        S.AmU      Filtered  state vector  A_t|t        (nobs+1 x m)
#        S.Pm       Predicted covariance of A_t|t-1      (nobs x m x m)
#        S.PmU      Filtered  covariance of A_t|t        (nobs+1 x m x m)
#        S.loglik   Value of likelihood function

# Output structure & dimensions
  n, m = size(Z);
  nobs  = size(Y,2);

  S.Am  = fill(NaN,m,nobs);   S.Pm  = fill(NaN,m,m,nobs);
  S.AmU = fill(NaN,m,nobs+1);   S.PmU = fill(NaN,m,m,nobs+1);
  S.loglik = 0;

  #______________________________________________________________________
  Au = A_0;  # A_0|0;
  Pu = P_0;  # P_0|0

  S.AmU[:,1]    = Au;
  S.PmU[:,:,1]  = Pu;

  for t = 1:nobs  #  iterates over time
      # A = A_t|t-1   & P = P_t|t-1
      A   = T*Au;
      P   = T*Pu*T' + Q;
      P   =  0.5 * (P+P');

      # handling the missing data
      y_t, Z_t, R_t, L_t = MissData(Y[:,t],Z,R);

      if isempty(y_t)
          Au = A;
          Pu = P;

      else
          PZ  = P*Z_t';
          iF  = inv(Z_t*PZ + R_t);
          PZF = PZ*iF;

          V   = y_t - Z_t*A;
          Au  = A  + PZF * V;
          Pu  = P  - PZF * PZ';
          Pu   =  0.5 * (Pu+Pu');
          S.loglik = S.loglik + 0.5*(log(det(iF))  - V'*iF*V);
      end

      S.Am[:,t]   = A;
      S.Pm[:,:,t] = P;

      # Au = A_t|t   & Pu = P_t|t

      S.AmU[:,t+1]    = Au;  # not re-writing Int64 for a variable declaration, but instead modifying elements of the matrix which i think is okay for im/mutability
      S.PmU[:,:,t+1]  = Pu;
  end # t
  if isempty(y_t)
      S.KZ = zeros(m,m);
  else
      S.KZ = PZF*Z_t;
  end

  return S
end

#______________________________________________________________________
function FIS!(S,Y,Z,R,T,Q)
#______________________________________________________________________
# Fixed intervall smoother (see Harvey, 1989, p. 154)
# FIS returns the smoothed state vector AmT and its covar matrix PmT
# Use this in conjnuction with function SKF
#______________________________________________________________________
# INPUT
#        Y         Data                                 (nobs x n)
#        S Estimates from Kalman filter SKF
#          S.Am   : Estimates     a_t|t-1                  (nobs x m)
#          S.Pm   : P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
#          S.AmU  : Estimates     a_t|t                    (nobs x m)
#          S.PmU  : P_t|t   = Cov(a_t|t)               (nobs x m x m)
# OUTPUT
#        S Smoothed estimates added to above
#          S.AmT  : Estimates     a_t|T                    (nobs x m)
#          S.PmT :  P_t|T   = Cov(a_t|T)               (nobs x m x m)
#          S.PmT_1 : Cov(a_ta_t-1|T)
#        where m is the dim of state vector and t = 1 ...T is time

  m, nobs        = size(S.Am);
  S.AmT           = zeros(m,nobs+1);
  S.PmT           = zeros(m,m,nobs+1);
  S.AmT[:,nobs+1]   = S.AmU[:,nobs+1]  ;
  S.PmT[:,:,nobs+1] = squeeze(S.PmU[:,:,nobs+1]);
  S.PmT_1[:,:,nobs] = (eye(m)-S.KZ) *T*squeeze(S.PmU[:,:,nobs]);

  J_2 = squeeze(S.PmU[:,:,nobs]) * T' * pinv(squeeze(S.Pm[:,:,nobs]));

  for t = nobs:-1:1
      PmU = squeeze(S.PmU[:,:,t]);
      Pm1 = squeeze(S.Pm[:,:,t]);
      P_T = squeeze(S.PmT[:,:,t+1]);
      P_T1 = squeeze(S.PmT_1[:,:,t]);

      J_1 = J_2;

      S.AmT[:,t]   = S.AmU[:,t] + J_1 * (S.AmT[:,t+1] - T * S.AmU[:,t]) ;
      S.PmT[:,:,t] = PmU        +  J_1 * (P_T - Pm1) * J_1';

     if t>1
          J_2 = squeeze(S.PmU[:,:,t-1]) * T' * pinv(squeeze(S.Pm[:,:,t-1]));
          S.PmT_1[:,:,t-1] = PmU*J_2'+J_1*(P_T1-T*PmU)*J_2';
      end
  end
  return S
end

#______________________________________________________________________
function MissData!(y,C,R)
#______________________________________________________________________
# PROC missdata
# PURPOSE: eliminates the rows in y & matrices Z, G that correspond to
#          missing data (NaN) in y
# INPUT    y             vector of observations at time t  (n x 1 )
#          S             KF system matrices             (structure)
#                        must contain Z & G
# OUTPUT   y             vector of observations (reduced)   (# x 1)
#          Z G           KF system matrices     (reduced)   (# x ?)
#          L             To restore standard dimensions     (n x #)
#                        where # is the nr of available data in y
#______________________________________________________________________
  ix = setdiff(find(isnan,y),1:length(y));
  e  = eye(size(y,1));  # identity matrix of size nxn == available data points in y
  L  = e[:,ix];         # identity truncated

  y  =  y[ix];  # useable data not equal to NaN
  C  =  C[ix,:];
  R  =  R[ix,ix]; # Revised variance matrix

  return y,C,R,L
end

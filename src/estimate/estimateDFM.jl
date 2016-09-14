#  Replication files for:
#  ""Nowcasting", 2010, (by Marta Banbura, Domenico Giannone and Lucrezia Reichlin),
# in Michael P. Clements and David F. Hendry, editors, Oxford Handbook on Economic Forecasting.
#
# The software can be freely used in applications.
# Users are kindly requested to add acknowledgements to published work and
# to cite the above reference in any resulting publications

# Make into subtype of AbstractFloat --> As in the DSGE Kalman code
type ResultsType{T} # NOT of the abstract type Model... X !<: AbstractModel{T}
    X_sm::Array{Float64,2}       #
    F::Array{Float64,2}       #
    C::Array{Float64,2}       #
    R::Array{Float64,2}       #
    A::Array{Float64,2}       #
    Q::Array{Float64,2}       #
    Mx::Array{Float64,2}      #
    Wx::Array{Float64,2}      #
    Z_0::Array{Float64,2}     #
    V_0::Array{Float64,2}     #
    r::Int                    #
    p::Int                    #
    Loglik::Float64           # Value of likelihood function
end


# Can explicitly state the type of Res if Res were to be only an input in these codes
# Because we are going to keep cycling back to Res, I want to break it out into its own type that we can use outside this function declaration
function estimateDFM(X::DataFrame,m::Model990)  # hand function Res if one exists and overwrite --> In main file will either load Res or initialize an empty Res
#                                      Res::DfmType  this is optional input
T,N = size(X)

r = m.r
p = m.p
nQ = length(find(m.nfrequency.=="q")) # calculate this value and store in the results vector, no need to supply this

i_idio = [trues(N-nQ,1);falses(nQ,1)] # mark the monthly variables with true logical

# blocks = Res.blocks
# blocks = m.dataBasket # symbols with integer array following (need to printout for each symbol using...)
blocks = falses(N,length(m.blockNames))
# initialize empty matrix and fill in the order of the X::DataFrame
idx_iQ = falses(N,1)
idx_iM = falses(N,1)
idx_iW = falses(N,1)
idx_iD = falses(N,1)
for (key, value) in m.dataBasket
    blocks[X.colindex[key],:] = value
    if m.nfrequency[key]=="q"
        idx_iQ[X.colindex[key],1] = true
    elseif m.nfrequency[key]=="m"
        idx_iM[X.colindex[key],1] = true
    elseif m.nfrequency[key]=="w"
        idx_iW[X.colindex[key],1] = true
    elseif m.nfrequency[key]=="d"
        idx_iD[X.colindex[key],1] = true
    end
end


#R*Lambda = q; Contraints on the loadings of the quartrly variables

R_mat = [  2 -1  0  0  0;
           3  0 -1  0  0;
           2  0  0 -1  0;
           1  0  0  0 -1]

q = zeros(4,1)

thresh = 1e-4

if nargin<4
    max_iter = 1000
end;
#--------------------------------------------------------------------------
# Preparation of the data
#--------------------------------------------------------------------------

#Removing missing values (for initial estimators)
#remove the leading and ending nans for the estimation
optionMethod = 2 # Remove leading and closing zeros
optionK = 3

if nargin<3
    # Standardise x
    Mx = mean(X)
    Wx = (std(X))
    xNaN = (X-repmat(Mx,T,1))./repmat(Wx,T,1)

    #--------------------------------------------------------------------------
    # Initial Conditions
    #--------------------------------------------------------------------------
    A, C, Q, R, Z_0, V_0 = InitCond(xNaN,r,p,blocks,optionMethod,optionK, R_mat,q,nQ,i_idio,X)

else

    A  = Res_old.A
    C  = Res_old.C
    Q  = Res_old.Q
    R  = Res_old.R
    Z_0 = Res_old.Z_0
    V_0 = Res_old.V_0
    Mx = Res_old.Mx
    Wx = Res_old.Wx
    # # Standardise x - if range of realtime is sufficiently less than range estimation, no need
    # Mx = mean(X);
    # Wx = (std(X));
    xNaN = (X-repmat(Mx,T,1))./repmat(Wx,T,1)

end;

# some auxiliary variables for the iterations
previous_loglik = -inf
num_iter = 0
LL = -inf
converged = 0

# y for the estimation is WITH missing data
y = xNaN'


#--------------------------------------------------------------------------
#THE EM LOOP
#--------------------------------------------------------------------------

#The model can be written as
#y = C*Z + e;
#Z = A*Z(-1) + v
#where y is NxT, Z is (pr)xT, etc


y_est = remNaNs_spline!(xNaN,optionMethod, optionK)';
LL = Array{Float64,1}()
DD = Array{Int,1}()
while (num_iter < max_iter) & ~converged
    C_new, R_new, A_new, Q_new, Z_0, V_0, loglik = EMstep(y_est, A, C, Q, R, Z_0, V_0, r,p,R_mat,q,nQ,i_idio,blocks);

    C = C_new
    R = R_new
    A = A_new
    Q = Q_new

    # Checking convergence
    if num_iter>2
        converged, decrease = em_converged(loglik, previous_loglik, thresh,1);
    end

    DD = [DD decrease]
    LL = [LL loglik]
    previous_loglik = loglik
    num_iter =  num_iter + 1

#     if mod(num_iter,10)==0
#         disp([num_iter loglik/(T*N)]);
#     end

end

#final run of the Kalman filter
Zsmooth = runKF(y, A, C, Q, R, Z_0, V_0)'
x_sm = Zsmooth[2:end,:]*C'

#--------------------------------------------------------------------------
#   Loading the structure with the results
#--------------------------------------------------------------------------
Res.X_sm = repmat(Wx,T,1).*x_sm+repmat(Mx,T,1)
Res.F = Zsmooth[2:end,:]
Res.C = C
Res.R = R
Res.A = A
Res.Q = Q
Res.Mx = Mx
Res.Wx = Wx
Res.Z_0 = Z_0
Res.V_0 = V_0
Res.r = r
Res.p = p
Res.Loglik = LL
  return Res
end

#--------------------------------------------------------------------------
#PROCEDURES
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
function em_converged!(loglik, previous_loglik, threshold, check_increased)
# EM_CONVERGED Has EM converged?
# [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
#
# We have converged if the slope of the log-likelihood function falls below 'threshold',
# i.e., |f(t) - f(t-1)| / avg < threshold,
# where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
# 'threshold' defaults to 1e-4.
#
# This stopping criterion is from Numerical Recipes in C p423
#
# If we are doing MAP estimation (using priors), the likelihood can decrase,
# even though the mode of the posterior is increasing.

if nargin < 3; threshold = 1e-4 end
if nargin < 4; check_increased = 1 end

converged = 0
decrease = 0

if check_increased
    if loglik - previous_loglik < -1e-3 # allow for a little imprecision
        print("*********likelihood decreased from $(previous_loglik) to $(loglik) \n")
        decrease = 1
    end
end

delta_loglik = abs(loglik - previous_loglik)
# avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2  # what is this eps value???
avg_loglik = (abs(loglik) + abs(previous_loglik))/2
if (delta_loglik / avg_loglik) < threshold; converged = 1 end

return converged, decrease
end

#--------------------------------------------------------------------------

function InitCond(x,r,p,blocks,optionMethod, optionK,Rcon,q,NQ,i_idio,idx_iM,idx_iQ)

pC = size(Rcon,2)
ppC = max(p,pC)
n_b = size(blocks,2)

xBal, indNaN = remNaNs_spline!(x,optionMethod,optionK)
T,N = size(xBal)
NM = N-NQ

xNaN = xBal
xNaN[indNaN] = NaN
C = []
A = []
Q = []
initV = []

res = xBal
resNaN = xNaN
indNaN[1:pC-1,:] = true
for i = 1:n_b
    r_i = r[i]
    #--------------------------------------------------------------------------
    # Observation equation
    #--------------------------------------------------------------------------
    C_i = zeros(N,r_i*ppC)
    idx_i = find(blocks[:,i])
    idx_iM = idx_i[idx_i<NM+1]
    idx_iQ = idx_i[idx_i>NM]

# TO be implemented when running a weekly and daily model
# Robust way to find frequency of the data for organizing prior to EMsteps
# idx_iMb = setdiff(find(idx_iM),idx_i)
# idx_iQb = setdiff(find(idx_iQ),idx_i)



    v, d = eigs(cov(res[:,idx_iMb]),r_i) # matlab options 'lm',OPTS
    C_i[idx_iMb,1:r_i] = v
    f = res[:,idx_iMb]*v
    F = []
    for kk = 0:max(p+1,pC)-1
        F = [F f[pC-kk:end-kk,:]]
    end
    Rcon_i = kron(Rcon,eye(r_i))
    q_i = kron(q,zeros(r_i,1))
    ff = F[:,1:r_i*pC]
    for j = idx_iQb'
        xx_j = resNaN[pC:end,j]
        if sum(~isnan(xx_j)) < size(ff,2)+2
            xx_j = res[pC:end,j]
        end
        ff_j = ff[~isnan(xx_j),:]
        xx_j = xx_j[~isnan(xx_j)]
        iff_j = inv(ff_j'*ff_j)
        Cc = iff_j*ff_j'*xx_j
        Cc = Cc - iff_j*Rcon_i'*inv(Rcon_i*iff_j*Rcon_i')*(Rcon_i*Cc-q_i)
        C_i[j,1:pC*r_i]=Cc'
    end
    ff = [zeros(pC-1,pC*r_i);ff]
    res = res - ff*C_i'
    resNaN = res
    resNaN(indNaN) = NaN
    C = [C C_i]

    #--------------------------------------------------------------------------
    # Transition equation
    #--------------------------------------------------------------------------
    z = F[:,1:r_i]
    Z = F[:,r_i+1:r_i*(p+1)]
    A_i = zeros(r_i*ppC,r_i*ppC)'
    A_temp = inv(Z'*Z)*Z'*z
    A_i[1:r_i,1:r_i*p] = A_temp'
    A_i[r_i+1:end,1:r_i*(ppC-1)] = eye(r_i*(ppC-1))

    ############
    Q_i = zeros(ppC*r_i,ppC*r_i)
    e = z  - Z*A_temp         # VAR residuals
    Q_i[1:r_i,1:r_i] = cov(e) # VAR covariance matrix

    initV_i = reshape(inv(eye((r_i*ppC)^2)-kron(A_i,A_i))*Q_i[:],r_i*ppC,r_i*ppC)

# Blk Diagonal only works for sparse matrices...?
    A = blkdiag(A,A_i)
    Q = blkdiag(Q,Q_i)
    initV = blkdiag(initV,initV_i)
end


R = diag(var(resNaN))

eyeN = eye(N)
# eyeN[:,~i_idio] = [];
eyeN = eyeN[:,i_idio]
# Initial conditions
C=[C eyeN]

ii_idio = find(i_idio)
n_idio = length(ii_idio)
B = zeros(n_idio)
S = zeros(n_idio)

for i = 1:n_idio
    R[ii_idio[i],ii_idio[i]] = 1e-04

    res_i = resNaN[:,ii_idio[i]]
    # number of leading zeros
    leadZero = max(find( [1:T]' .== cumsum(isnan(res_i))))
    endZero = max(find( [1:T]' .== cumsum(isnan(res_i[end:-1:1]))))

    res_i = res[:,ii_idio[i]]
    # res_i[end-endZero:endZero] = [];
    # res_i[end-endZero:end] = [];  % Was this an error in original code ==> the prior line code makes no sense (eliminate 1 row at endZero) --> works fine if endZero==0
    res_i = res_i[leadZero+1:end-endZero] # endZero-1
    # res_i[1:leadZero] = [];

    BM[i,i] = inv(res_i[1:end-1]'*res_i[1:end-1])*res_i[1:end-1]'*res_i[2:end,:]
    SM[i,i] = cov(res_i[2:end]-res_i[1:end-1]*B[i,i])
end
initViM = diag(1./diag(eye(size(BM,1))-BM.^2)).*SM


C = [C [zeros(NM,5*NQ);kron(eye(NQ),[1 2 3 2 1])]]
Rdiag = diag(R)
sig_e = Rdiag(NM+1:N)/19
Rdiag[NM+1:N] = 1e-04
R = diag(Rdiag)


rho0 = 0.1
BQ = kron(eye(NQ),[[rho0 zeros(1,4)];[eye(4),zeros(4,1)]])
temp = zeros(5)
temp[1,1] = 1
SQ = kron(diag((1-rho0^2)*sig_e),temp)

initViQ = reshape(inv(eye((5*NQ)^2)-kron(BQ,BQ))*SQ[:],5*NQ,5*NQ)

A = blkdiag(A, BM, BQ)
Q = blkdiag(Q, SM, SQ)

# Initial conditions
initZ = zeros(size(A,1),1)
initV = blkdiag(initV, initViM, initViQ)

    return A, C, Q, R, initZ, initV
end


#--------------------------------------------------------------------------
function EMstep(y, A, C, Q, R, Z_0, V_0, r,p,R_mat,q,nQ,i_idio,blocks)

n,T = size(y)
nM = n-nQ

pC = size(R_mat,2)
ppC = max(p,pC)


n_b = size(blocks,2)

# Compute the (expected) sufficient statistics for a single Kalman filter sequence.

#Running the Kalman filter with the current estimates of the parameters
Zsmooth, Vsmooth, VVsmooth, loglik = runKF(y, A, C, Q, R, Z_0, V_0)

A_new = A
Q_new = Q
V_0_new = V_0

for i = 1:n_b
    r_i = r[i]
    rp = r_i*p
    rp1 = sum(r[1:i-1])*ppC

    A_i = A[rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC]
    Q_i = Q[rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC]

    EZZ = Zsmooth[rp1+1:rp1+rp,2:end]*Zsmooth[rp1+1:rp1+rp,2:end]'
        +sum(Vsmooth[rp1+1:rp1+rp,rp1+1:rp1+rp,2:end],3)                       #E(Z'Z)
    EZZ_BB = Zsmooth[rp1+1:rp1+rp,1:end-1]*Zsmooth[rp1+1:rp1+rp,1:end-1]'
        +sum(Vsmooth[rp1+1:rp1+rp,rp1+1:rp1+rp,1:end-1],3) #E(Z(-1)'Z_(-1))
    EZZ_FB = Zsmooth[rp1+1:rp1+rp,2:end]*Zsmooth[rp1+1:rp1+rp,1:end-1]'
        +sum(VVsmooth[rp1+1:rp1+rp,rp1+1:rp1+rp,:],3) #E(Z'Z_(-1))

    A_i[1:r_i,1:rp] = EZZ_FB[1:r_i,1:rp] * inv(EZZ_BB[1:rp,1:rp])
    Q_i[1:r_i,1:r_i] = (EZZ[1:r_i,1:r_i] - A_i[1:r_i,1:rp]*EZZ_FB[1:r_i,1:rp]') / T

    A_new[rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC] = A_i
    Q_new[rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC] = Q_i
    V_0_new[rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC] = Vsmooth[rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC,1];
end

rp1 = sum(r)*ppC;
niM = sum(i_idio[1:nM]);
# idiosyncratic
EZZ = diag(diag(Zsmooth[rp1+1:end,2:end]*Zsmooth[rp1+1:end,2:end]'))
    +diag(diag(sum(Vsmooth[rp1+1:end,rp1+1:end,2:end],3)))                        #E(Z'Z)
EZZ_BB = diag(diag(Zsmooth[rp1+1:end,1:end-1]*Zsmooth[rp1+1:end,1:end-1]'))
    +diag(diag(sum(Vsmooth[rp1+1:end,rp1+1:end,1:end-1],3)))  #E(Z(-1)'Z_(-1))
EZZ_FB = diag(diag(Zsmooth[rp1+1:end,2:end]*Zsmooth[rp1+1:end,1:end-1]'))
    +diag(diag(sum(VVsmooth[rp1+1:end,rp1+1:end,:],3)))  #E(Z'Z_(-1))

A_i = EZZ_FB * diag(1./diag((EZZ_BB)))
Q_i = (EZZ - A_i*EZZ_FB') / T

A_new[rp1+1:rp1+niM,rp1+1:rp1+niM] = A_i[1:niM,1:niM]
Q_new[rp1+1:rp1+niM,rp1+1:rp1+niM] = Q_i[1:niM,1:niM]
V_0_new[rp1+1:rp1+niM,rp1+1:rp1+niM] = diag(diag(Vsmooth[rp1+1:rp1+niM,rp1+1:rp1+niM,1]))

Z_0 = Zsmooth[:,1]  #zeros(size(Zsmooth,1),1); #


nanY = isnan(y)
y[nanY] = 0

# LOADINGS
C_new = C


# Not really understanding the sparse matrix getting generated below, so trying to make 1-1 with matlab syntax
############
# Blocks
bl = unique(blocks, 1)
n_bl = size(bl,1);
bl_idxM = Array{Int,2}() # could try initializing to Boolto speed up code
bl_idxQ = Array{Int,2}()
# bl_idxM = []
# bl_idxQ = []
R_con = []
q_con = []
for i = 1:n_b
    bl_idxQ = [bl_idxQ repmat(bl[:,i],1,r[i]*ppC)]
    bl_idxM = [bl_idxM repmat(bl[:,i],1,r[i]) zeros(n_bl,r[i]*(ppC-1))]
    R_con = blkdiag(R_con, kron(R_mat,eye(r[i])))
    q_con = [q_con; zeros(r[i]*size(R_mat,1),1)]
end

# Don't think we need as Integers 0,1 are treated as boolean anyway
# bl_idxM = logical(bl_idxM);
# bl_idxQ = logical(bl_idxQ);

#idio
i_idio_M = i_idio[1:nM]
n_idio_M = length(find(i_idio_M))
c_i_idio = cumsum(i_idio)

for i = 1:n_bl
    bl_i = bl[i,:]
    rs = sum(r[bl_i])

    # idx_i = find(ismember(blocks,bl_i,'rows'))
    # Not sure Julia equivalent
    idx_i = []
    for iii_ = 1:size(blocks,1)
        if isequal(blocks[iii_,:],bl_i)
            idx_i = [idx_i; iii_]
        end
    end


    # MONTHLY
    idx_iM = idx_i[idx_i<nM+1]
    n_i = length(idx_iM)

    denom = zeros(n_i*rs,n_i*rs)
    nom = zeros(n_i,rs)

    i_idio_i = i_idio_M(idx_iM)
    i_idio_ii = c_i_idio(idx_iM)
    i_idio_ii = i_idio_ii(i_idio_i)
    for t=1:T
        nanYt = diag(~nanY[idx_iM,t])
        denom = denom + kron(Zsmooth[bl_idxM[i,:],t+1]*Zsmooth[bl_idxM[i,:],t+1]'
            +Vsmooth[bl_idxM[i,:],bl_idxM[i,:],t+1],nanYt)
        nom = nom + y[idx_iM,t]*Zsmooth[bl_idxM[i,:],t+1]' #here's the modification
            -nanYt[:,i_idio_i]*(Zsmooth[rp1+i_idio_ii,t+1]*Zsmooth[bl_idxM[i,:],t+1]'
            +Vsmooth[rp1+i_idio_ii,bl_idxM[i,:],t+1])
    end
    vec_C = inv(denom)*nom[:]
    C_new[idx_iM,bl_idxM[i,:]] = reshape(vec_C,n_i,rs)

    # QUARTERLY
   idx_iQ = idx_i[idx_i>nM]
   rps = rs*ppC

   R_con_i = R_con[:,bl_idxQ[i,:]]
   q_con_i = q_con
   no_c = ~(any(R_con_i,2))
   R_con_i[no_c,:] = []
   q_con_i[no_c,:] = []


   for j = idx_iQ'
       denom = zeros(rps,rps)
       nom = zeros(1,rps)
       idx_jQ = j-nM
       i_idio_jQ = [rp1+n_idio_M+5*(idx_jQ-1)+1:rp1+n_idio_M+5*idx_jQ]
       V_0_new[i_idio_jQ,i_idio_jQ] = Vsmooth[i_idio_jQ,i_idio_jQ,1]
       A_new[i_idio_jQ[1],i_idio_jQ[1]] = A_i[i_idio_jQ[1]-rp1,i_idio_jQ[1]-rp1]
       Q_new[i_idio_jQ[1],i_idio_jQ[1]] = Q_i[i_idio_jQ[1]-rp1,i_idio_jQ[1]-rp1]

       for t=1:T
           nanYt = diag(~nanY[j,t])
           denom = denom + kron(Zsmooth[bl_idxQ[i,:],t+1]*Zsmooth[bl_idxQ[i,:],t+1]'
                +Vsmooth[bl_idxQ[i,:],bl_idxQ[i,:],t+1],nanYt)
            nom = nom + y[j,t]*Zsmooth[bl_idxQ[i,:],t+1]'
            nom = nom -
                nanYt*([1 2 3 2 1]*Zsmooth[i_idio_jQ,t+1]*Zsmooth[bl_idxQ[i,:],t+1]'+
                [1 2 3 2 1]*Vsmooth[i_idio_jQ,bl_idxQ[i,:],t+1])
        end

        C_i = inv(denom)*nom'
        C_i_constr = C_i - inv(denom)*R_con_i'*inv(R_con_i*inv(denom)*R_con_i')*(R_con_i*C_i-q_con_i)
        C_new[j,bl_idxQ[i,:]] = C_i_constr

    end
end


R_new = zeros(n,n)
for t=1:T
    nanYt = diag(~nanY[:,t])
    R_new = R_new + (y[:,t]-nanYt*C_new*Zsmooth[:,t+1])*(y[:,t]-nanYt*C_new*Zsmooth[:,t+1])'
        +nanYt*C_new*Vsmooth[:,:,t+1]*C_new'*nanYt
        +(eye(n)-nanYt)*R*(eye(n)-nanYt)
end

R_new = R_new/T
RR = diag(R_new)  #RR(RR<1e-2) = 1e-2;
RR[i_idio_M] = 1e-04
RR[nM+1:end] = 1e-04
R_new = diag(RR)

return C_new, R_new, A_new, Q_new, Z_0, V_0, loglik
end

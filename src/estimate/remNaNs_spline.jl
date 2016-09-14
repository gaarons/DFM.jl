#  Replication files for:
#  ""Nowcasting", 2010, (by Marta Banbura, Domenico Giannone and Lucrezia Reichlin),
# in Michael P. Clements and David F. Hendry, editors, Oxford Handbook on Economic Forecasting.
#
# The software can be freely used in applications.
# Users are kindly requested to add acknowledgements to published work and
# to cite the above reference in any resulting publications
function remNaNs_spline!(X::Array{Float64,2},optionMethod::Int, k::Int) # overwrite X

# put outside with the module declaration
# using Interpolations

T,N=size(X);
# create an indNaN map for the rest of the function
# indNaN = fill(false,T,N)
# indNaN[find(isnan,X)] = true
indNaN = isnan(X)
# Changed only option 1 and 2, re-created the functions nanmedian, nanmean
#TODO change the rest of the options
if optionMethod == 1
        for i = 1:N
            x = X[:,i]
            x[indNaN[:,i]]=nanmedian(x)

            # User defined MA filter
            x_preMA = [mean(x[i-k:i+k]) for i = k+1:length(x)-k]
            x_MA[k+1:end] = x_preMA

            # x_MA =filter(ones(2k+1,1)/(2k+1), 1, [x(1)*ones(k,1);x;x[end]*ones(k,1)])
            # x_MA=x_MA[2k+1:end]

            x[indNaN[:,i]]=x_MA[indNaN[:,i]]
            X[:,i]=x
        end
elseif optionMethod == 2 # replace missing values after removing leading and closing zeros
        rem1=(sum(indNaN,2)>N*0.8)  # observations where 80% of the set is NaN [nx1] = [false/true]

        nanLead=(cumsum(rem1,1).==[1:T])
        nanEnd=(cumsum(rem1[end:-1:1]).==[1:T])
        nan_temp=nanEnd[end:-1:1]
        nanEnd=nan_temp; # extra step b/c of julia
        nanLE=(nanLead | nanEnd)
        X_temp = X[find(~nanLE),:]
        X = X_temp
        # indNaN = fill(false,T,N)
        # indNaN[find(isnan,X)] = true
        indNaN = isnan(X)
        for i = 1:N
            x = X[:,i]
        #     isnanx = fill(false,T,N)
        #     isnanx[find(isnan,X)] = true
            indNaN = isnan(X)
            t1 = min(find(~isnanx))
            t2 = max(find(~isnanx))
            # not-a-knot ensure the y''' = 0 = y" = y'
            # x[t1:t2] = spline(find(~isnanx),x[~isnanx],[t1:t2]);
            # x[t1:t2] = spline(find(~isnanx),x[~isnanx],[t1:t2])

            # from interpolate package -- # This algorithm ensures that = 0 = y" = y' but y''' not set to zero
            # no functionality for extrapolating
            # x[t1:t2] = interpolate(x[~isnanx], BSpline(Cubic(Natural())), OnGrid())

# Dierckk Package for interpolations with cubic spline and not-a-knot derivatives = 0
spl_obj = Spline1D(find(~isnanx), x[~isnanx], find(~isnanx), w=ones(length(x)), k=3, bc="nearest")
x[t1:t2] = evaluate(spl_obj,[t1:t2])
        #     isnanx = fill(false,T,N)
        #     isnanx[find(isnan,X)] = true

            # Think this was incorrect in old code --> # indNaN = isnan(X)
            indNaN = isnan(x)
            x[isnanx]=median(x) # outside the leading and ending NaN use median
            # Simple moving average filter for the median set at the front and back of the time-series
            # Window-size set to 2k


            x_preMA = [mean(x[i-k:i+k]) for i = k+1:length(x)-k]
            x_MA[k+1:end] = x_preMA
            # x_MA =filter (ones(2k+1,1)/(2k+1),   1,   [x[1]*ones(k,1);x;x[end]*ones(k,1)]);
            # x_MA=x_MA[2k+1:end];
            x[isnanx]=x_MA[isnanx]  # Still using the median(x) values at the begining of the x vector...?
            X[:,i]=x
        end
    # elseif optionMethod == 3 # only remove rows with leading and closing zeros
    #     rem1=(sum(indNaN,2)==N);
    #     nanLead=(cumsum(rem1)==(1:T)');
    #     nanEnd=(cumsum(rem1(end:-1:1))==(1:T)');
    #     nanEnd=nanEnd(end:-1:1);
    #     nanLE=(nanLead | nanEnd);
    #     X(nanLE,:)=[];
    #     indNaN = isnan(X);
    # elseif optionMethod == 4 # remove rows with leading and closing zeros & replace missing values
    #
    #     rem1=(sum(indNaN,2)==N);
    #     nanLead=(cumsum(rem1)==(1:T)');
    #     nanEnd=(cumsum(rem1(end:-1:1))==(1:T)');
    #     nanEnd=nanEnd(end:-1:1);
    #     nanLE=(nanLead | nanEnd);
    #     X(nanLE,:)=[];
    #     indNaN=isnan(X);
    #     for i = 1:N
    #         x = X[:,i);
    #         isnanx = isnan(x);
    #         t1 = min(find(~isnanx));
    #         t2 = max(find(~isnanx));
    #         x(t1:t2) = spline(find(~isnanx),x(~isnanx),(t1:t2)');
    #         isnanx = isnan(x);
    #         x(isnanx)=nanmedian(x);
    #         x_MA =filter (ones(2*k+1,1)/(2*k+1),1,[x(1)*ones(k,1);x;x(end)*ones(k,1)]);
    #         x_MA=x_MA(2*k+1:end);
    #         x(isnanx)=x_MA(isnanx);
    #         X[:,i)=x;
    #     end
    # elseif optionMethod == 5 # replace missing values
    #     indNaN=isnan(X);
    #     for i = 1:N
    #         x = X[:,i);
    #         isnanx = isnan(x);
    #         t1 = min(find(~isnanx));
    #         t2 = max(find(~isnanx));
    #         x(t1:t2) = spline(find(~isnanx),x(~isnanx),(t1:t2)');
    #         isnanx = isnan(x);
    #         x(isnanx)=nanmedian(x);
    #         x_MA =filter (ones(2*k+1,1)/(2*k+1),1,[x(1)*ones(k,1);x;x(end)*ones(k,1)]);
    #         x_MA=x_MA(2*k+1:end);
    #         x(isnanx]=x_MA(isnanx);
    #         X[:,i]=x;
        end
return X, indNaN
end

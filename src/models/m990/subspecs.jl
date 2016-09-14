"""
`init_subspec!(m::Model990)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::Model990)
    # Observable::Symbol,
    # transformation::String
    # nfrequency::String
    # startDate::String
    # endDate::String
    # Factor Block Structure       [Global, Macro, Financial]
    if subspec(m) == "ss0"
        # ss0(m::Model990)
        m.p = 1
        m.r = [1 1 1]
        m.blockNames = ["Global" "Macro" "Financial"]
        observables =  [:GDP          "1" "m" "" "" 1 1 0; # quarterly output growth
                        :GDPCTPI      "1" "m" "" "" 1 0 0;
                        :PCE          "1" "m" "" "" 1 1 0;
                        :FPI          "1" "m" "" "" 1 0 1;
                        :obs_wages    "1" "m" "" "" 1 1 0]
    elseif subspec(m) == "ss2"
        # ss0(m::Model990)`
        m.p = 1
        m.r = [1 1 1]
        m.blockNames = ["Global" "Macro" "Financial"]
        observables =  [:GDP          "1" "m" "" "" 1 1 0; # quarterly output growth
                        :GDPCTPI      "1" "m" "" "" 1 0 0;
                        :PCE          "1" "m" "" "" 1 1 0;
                        :FPI          "1" "m" "" "" 1 0 1;
                        :obs_wages    "0" "m" "" "" 1 1 0]
    else
        error("This subspec is not defined.")
    end

    # Breakout into individual maps --> we keep subspec format as in MATLAB b/c easier to read
    for (i,k) in enumerate(observables[:,1]);  m.observables[k] = i end
    for (i,k) in enumerate(observables[:,1]);  m.transf[k] = observables[i,2] end
    for (i,k) in enumerate(observables[:,1]);  m.nfrequency[k] = observables[i,3] end
    for (i,k) in enumerate(observables[:,1]);  m.startDate[k] = observables[i,4] end
    for (i,k) in enumerate(observables[:,1]);  m.endDate[k] = observables[i,5] end
    for (i,k) in enumerate(observables[:,1]);  m.dataBasket[k] = observables[i,6:end] end
    # numBaskets = size(observables,2) - 5;  # Eventually want Dict creation and basketing dynamic w.r.t. however many baskets are called.
    m <= Setting(:num_baskets, size(observables,2) - 5, "Number of factor baskets")

# determined elsewhere    m <= Setting(:frequency, "m", "Frequency of the model")

# I want to take the approach of collecting symbols in each block type if the boolean is accepted
# for ii_ = 1:1 #numBaskets
#     # for (i,k) in enumerate(observables[:,1])
#     #       if (observables[i,6+ii_]==1) m.("block_$(ii_)")[k] = (observables[i,6+ii_]==1)
#     #       end
#     # end
#
#     # for (i,k) in enumerate(observables[:,1]);  m.gBasket[k] = (observables[i,6]==1) end
#     # for (i,k) in enumerate(observables[:,1]);  m.mBasket[k] = (observables[i,7]==1) end
#     # for (i,k) in enumerate(observables[:,1]);  m.fBasket[k] = (observables[i,8]==1) end
# end

end # function

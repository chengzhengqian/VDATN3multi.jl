# collection of field to define the model

mutable struct Model
    # N=2, or N=3 (it is trivial to add N=1,HF)
    N_time_step::Integer
    # number the spin orbitals
    N_spin_orbital::Integer
    # symmetry of the spin orbitals
    symmetry
    # if we want to constrol the density, determines the density, emtpy indicates that there is no constraint on density
    n_target
    # local variational parameters, if n_target is emtpy, then w_para also determined the density
    w_para
    # emtpy for N=2
    G12_para
    # N_symmetry or N_symmtry*2, (in the former case, β< = β>; in the latter case, β<,β>), emtpy for N=2
    b_para
    # interaction, two particle [(idx1,idx2,c)...]
    interaction
    # chemcial potential, to control the density, [(idx1,c)...], empty for fixed density case
    chemical_potential
    # lattice part, N_symmetry, [e_fn...], for each e_fn (for a given symmetry group,), e_fn(p)  gives the energy for p in [0,1], where p is the accumulated probability
    e_fns
    # options, cutoff and N_k_points ect
    options::Dict
    # results
    obs::Dict
end

# we can define macro for this process
# "1"*"2"
macro define_option(name)
    # we can use mutiple dispatch. Also consider add type for the setter in future
    fn_get=Symbol("option_"*name)
    fn_set=Symbol("option_"*name)
    Model=Symbol("Model")
    quote
        function $(esc(fn_get))(model::$(esc(Model)))
            model.options[$(name)]
        end        
        function $(esc(fn_set))(model::$(esc(Model)),val)
            model.options[$(name)]=val
        end        
    end    
end

# """
# add variable to model.obs, using the variable's name as key
# String(:a)
# """
# macro set_obs(model,name)
#     quote
#         $(esc(model)).obs[$(String(name))]=$(esc(name))
#     end    
# end
# !! we have made the general version
# we may need need to set multiple names
# expr=quote
#     model.obs["y"]=y
#     model.obs["x"]=x
# end
# dump(expr)
macro set_obs(m,name...)
    exprs=[]
    for (idx_,name_) in enumerate(name)
        push!(exprs,LineNumberNode(idx_))
        push!(exprs,Expr(:escape,:($(m).obs[$(String(name_))]=$(name_))))
    end
    Expr(:block,exprs...)
end

macro get_obs(m,name...)
    exprs=[]
    for (idx_,name_) in enumerate(name)
        push!(exprs,LineNumberNode(idx_))
        push!(exprs,Expr(:escape,:($(name_)=$(m).obs[$(String(name_))])))
    end
    Expr(:block,exprs...)
end


# m=:model
# name=[:a,:b]
# dump(a)
# dump(b)
# dump(a.args[2])
# typeof(a.args[1])
# @macroexpand(@set_obs model a d c)
# @macroexpand(@get_obs model a d c)
# @macroexpand(@get_obs model a )
# @set_obs model a d
# @macroexpand(@set_obs model a)
# wether we need to compute eασ, true, for fixed the density case and first run

"""
customize how model looks in repl
"""
function Base.show(io::IO, model::Model)
    compact = get(io, :compact, false)
    if !compact
        print("<model N=$(model.N_time_step) N_spin_orb=$(model.N_spin_orbital) fixed_n=$(option_is_density_fixed(model))>")
    else
        print("<model N=$(model.N_time_step)>")
    end    
end


# ---
# title: Simulation de la végétation d’un corridor sous ligne électrique
# repository: tpoisot/BIO245-modele
# auteurs:
#    - nom: Mouchaimech
#      prenom: Kazem
#      matricule: 20232897
#      github: Kaz1711
#    - nom: Nivelle
#      prenom: Anna
#      matricule: 20349280
#      github: anna2808-prog
# ---

# # Introduction
# L'aménagement de corridors sous les lignes électriques nécessite un compromis
# entre la sécurité et des infrastructures et le maintien de la biodiversité.
# Une végétation trop haute représente un risque pour les lignes, alors qu'un 
# corridor trop dégagé réduit la biodiversité. 
# Dans ce travail, nous utilisons un modèle de transitions végétales pour simuler
# l'évolution d'un corridor composé de 200 parcelles. le modèle permet d'évaluer 
# si une intervention initiale permet d'atteindre les objectifs. 20% des parcelles
# doivent être végétalisées, dont 30% herbes et 70% buissons. Et que la variété de
# la moins abondante doit représenter au moins 30% des buissons. 
# # Présentation du modèle

# # Implémentation
# ## Packages nécessaires
using CairoMakie
using Distributions

import Random
Random.seed!(2045)
# ## Fonctions
function check_transition_matrix!(T)
    for ligne in axes(T, 1)
        if sum(T[ligne, :]) != 1
            @warn "La somme de la ligne $(ligne) n'est pas égale à 1 et a été modifiée"
            T[ligne, :] ./= sum(T[ligne, :])
        end
    end
    return T
end

function check_function_arguments(transitions, states)
    if size(transitions, 1) != size(transitions, 2)
        throw("La matrice de transition n'est pas carrée")
    end

    if size(transitions, 1) != length(states)
        throw("Le nombre d'états ne correspond psa à la matrice de transition")
    end
    return nothing
end

function _sim_stochastic!(timeseries, transitions, generation)
    for state in axes(timeseries, 1)
        pop_change = rand(Multinomial(timeseries[state, generation], transitions[state, :]))
        timeseries[:, generation+1] .+= pop_change
    end
end

function _sim_determ!(timeseries, transitions, generation)
    pop_change = (timeseries[:, generation]' * transitions)'
    timeseries[:, generation+1] .= pop_change
end

function simulation(transitions, states; generations=500, stochastic=false)

    check_transition_matrix!(transitions)
    check_function_arguments(transitions, states)

    _data_type = stochastic ? Int64 : Float32
    timeseries = zeros(_data_type, length(states), generations + 1)
    timeseries[:, 1] = states

    _sim_function! = stochastic ? _sim_stochastic! : _sim_determ!

    for generation in Base.OneTo(generations)
        _sim_function!(timeseries, transitions, generation)
    end

    return timeseries
end

# ## States
# ## Barren, Grass, Shrub_1, Shrub_2
s = [150, 20, 10, 10]
states = length(s)
patches = sum(s)

# ## Transitions
T = zeros(Float64, states, states)
T[1, :] = [0.97, 0.01, 0.01, 0.01] # vide reste souvent vide ou devient herbe
T[2, :] = [0.12, 0.82, 0.03, 0.03] # herbe peut devenir buisson
T[3, :] = [0.09, 0.02, 0.87, 0.02] # les buissons sont stables
T[4, :] = [0.09, 0.02, 0.02, 0.87] # les buissons sont stables
T

println("Somme de chaque ligne :")
println(sum(T, dims=2))

states_names = ["Barren", "Grasses", "Shrub_1", "Shrub_2"]
states_colors = [:grey40, :orange, :teal, :blue]

# ##Simulations et presentation des résultats

f = Figure()
ax = Axis(f[1, 1], xlabel="Nb. générations", ylabel="Nb. parcelles")

# ##Stochastic simulation
for _ in 1:100
    sto_sim = simulation(T, s; stochastic=true, generations=200)
    for i in eachindex(s)
        lines!(ax, sto_sim[i, :], color=states_colors[i], alpha=0.1)
    end
end

# ##Deterministic simulation
det_sim = simulation(T, s; stochastic=false, generations=200)
for i in eachindex(s)
    lines!(ax, det_sim[i, :], color=states_colors[i], alpha=1, label=states_names[i], linewidth=4)
end


axislegend(ax)
tightlimits!(ax)
current_figure()

# ## Verifications de l'équilibre
success = 0

for i in 1:100
    
    sim = simulation(T, s; stochastic=true, generations=200)
    final = sim[:, end]

    vegetation = final[2] + final[3] + final[4]
    shrubs = final[3] + final[4]

    cond1 = abs(vegetation-40) <= 5
    cond2 = abs(final[2] - 12) <= 5
    cond3 = min(final[3], final[4]) >= 0.3*shrubs

    if cond1 && cond2 && cond3
        success += 1
    end
end

println("Success rate = ", success/100)
# # Discussion

# On peut aussi citer des références dans le document `references.bib`,
# @ermentrout1993cellular -- la bibliographie sera ajoutée automatiquement à la
# fin du document.

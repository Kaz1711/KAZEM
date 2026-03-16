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
# Nous utilisons un modèle de transition décrivant l’évolution de l’état de parcelles de végétation au cours du temps. Chaque parcelle peut se trouver dans l’un 
# des quatre états suivants : sol nu (barren), herbes (grasses), buisson de type 1 (shrub_1) ou buisson de type 2 (shrub_2). Les changements d’état entre les 
# générations sont décrits par une matrice de transition, dont les lignes correspondent de la parcelle à chaque génération et les colonnes à la probabilité de 
# transition vers un autre état lors de la génération suivante. Les probabilités sont normalisées pour que la somme de chaque ligne soit égale à 1. Le modèle est 
# simulé à la fois de façon déterministe et de façon stochastique. Dans le cas déterministe, les transitions sont calculées directement à partir des probabilités 
# de la matrice. Dans le cas stochastique, les transitions sont tirées aléatoirement à l’aide d’une distribution multinomiale pour introduire une variabilité 
# entre les simulations. 

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
# ## ajouter un deuxième état de buisson pour avoior quatre états
# ## Barren, Grass, Shrub_1, Shrub_2
s = [150, 0, 25, 25] # un quatrième état de buisson est ajouté
states = length(s)
patches = sum(s)

# ## Transitions
# ## ajouter un quatrième ligne et colonne pour la matrice de transition pour le deuxième type de buisson
T = zeros(Float64, states, states)
T[1, :] = [0.97, 0.01, 0.01, 0.01] # vide reste souvent vide avec une petite chance de devenir herbe ou buisson
T[2, :] = [0.11, 0.84, 0.02, 0.03] # herbe reste souvent herbe avec une chance de devenir vide ou buisson
T[3, :] = [0.13, 0.01, 0.84, 0.02] # les buissons restent souvent des buissons avec une chance de devenir vide ou herbe ou se transformer en un autre type
T[4, :] = [0.13, 0.01, 0.02, 0.84] # les buissons restent souvent des buissons avec une chance de devenir vide ou herbe ou se transformer en un autre type 
T

println("Somme de chaque ligne :")
println(sum(T, dims=2))

states_names = ["Barren", "Grasses", "Shrub_1", "Shrub_2"] # ajouter le nom du deuxième type de buisson
states_colors = [:grey40, :orange, :teal, :blue] #ajouter une couleur pour le deuxième type de buisson

# ##Simulations et presentation des résultats

f = Figure()
ax = Axis(f[1, 1], xlabel="Nb. générations", ylabel="Nb. parcelles")

# ##Stochastic simulation
for _ in 1:100
    local sto_sim = simulation(T, s; stochastic=true, generations=200)
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
display(f) # pour afficher la figure

# ## Verifications de l'équilibre
function check_success(T, s)
    # on definit une fraction de succès pour les simulations stochastiques
    success = 0 
    for i in 1:100
        sim = simulation(T, s; stochastic=true, generations=200)
        final = sim[:, end] # on regarde la composition du corridor à la fin de la simulation
        vegetation = final[2] + final[3] + final[4] # on calcule la quantité de végétation totale
        shrubs = final[3] + final[4] # on calcule la quantité de buissons totale
        cond1 = abs(vegetation-40) <= 5 # on vérifie que la quantité de végétation est proche de 20% du corridor
        cond2 = abs(final[2] - 12) <= 5 # on vérifie que la quantité d'herbe est proche de 6% du corridor
        cond3 = min(final[3], final[4]) >= 0.3*shrubs # on vérifie que la variété de la moins abondante ne représente pas moins que 30% des buissons

        if cond1 && cond2 && cond3 # si les trois conditions sont vérifiées, on considère que la simulation est un succès
        success += 1 # on additionne 1 au nombre de succès
        end
     end
     return success/100 # on retourne la fraction de succès
end

println("Success rate = ", check_success(T, s)) # on montre la fraction de succès des simulations stochastiques
# # Discussion

# On peut aussi citer des références dans le document `references.bib`,
# @ermentrout1993cellular -- la bibliographie sera ajoutée automatiquement à la
# fin du document.

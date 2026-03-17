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
# des quatre états suivants : sol nu (Barren), herbes (Grasses), buisson de type 1 (Shrub_1) ou buisson de type 2 (Shrub_2). Les changements d’état entre les 
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
# vérifier que la somme des probabilités est égale à 1
function check_transition_matrix!(T)
    # les valeurs sont normalisées si la somme des probabilités n'est pas égale à 1
    for ligne in axes(T, 1)
        if sum(T[ligne, :]) != 1
            @warn "La somme de la ligne $(ligne) n'est pas égale à 1 et a été modifiée"
            T[ligne, :] ./= sum(T[ligne, :])
        end
    end
    return T
end

# vérifier que la matrice de transition est carrée et que le nombre d'états lui correspond
function check_function_arguments(transitions, states)
    if size(transitions, 1) != size(transitions, 2)
        throw("La matrice de transition n'est pas carrée")
    end

    if size(transitions, 1) != length(states)
        throw("Le nombre d'états ne correspond pas à la matrice de transition")
    end
    return nothing
end

# simulation stochastique
function _sim_stochastic!(timeseries, transitions, generation)
    for state in axes(timeseries, 1)
        # tirage aléatoire du nombre de transitions par état
        pop_change = rand(Multinomial(timeseries[state, generation], transitions[state, :]))
        timeseries[:, generation+1] .+= pop_change
    end
end

# simulation déterministe
function _sim_determ!(timeseries, transitions, generation)
    pop_change = (timeseries[:, generation]' * transitions)'
    timeseries[:, generation+1] .= pop_change
end

# réaliser les simulations en vérifiant les paramètres de la matrice de transition
function simulation(transitions, states; generations=500, stochastic=false)

    check_transition_matrix!(transitions)
    check_function_arguments(transitions, states)

    # les parcelles sont entières si modèle déterministe et fractionnées si stochastique
    _data_type = stochastic ? Int64 : Float32
    timeseries = zeros(_data_type, length(states), generations + 1)
    timeseries[:, 1] = states

    # choix de la simulation
    _sim_function! = stochastic ? _sim_stochastic! : _sim_determ!
    # calculer les générations suivantes
    for generation in Base.OneTo(generations)
        _sim_function!(timeseries, transitions, generation)
    end

    return timeseries
end

# ## States
# ## Barren, Grass, Shrub_1, Shrub_2
s = [175, 0, 8, 17] 
states = length(s)
patches = sum(s)

# ## Transitions
# ## ajouter un quatrième ligne et colonne pour la matrice de transition pour le deuxième type de buisson
T = zeros(Float64, states, states)
T[1, :] = [0.97, 0.01, 0.01, 0.01] # vide reste souvent vide 
T[2, :] = [0.11, 0.84, 0.02, 0.03] # herbe reste surtout herbe avec une chance de devenir vide 
T[3, :] = [0.13, 0.01, 0.84, 0.02] # buissons restent souvent des buissons avec une chance de devenir vide 
T[4, :] = [0.13, 0.01, 0.02, 0.84]  
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
display(f)
save("travail-1.png", f) # sauvegarder la figure pour une utilisation dans le rapport 

# ## Verifications de l'équilibre
function check_success(T, s)
    # obtenir une fraction de succès pour les simulations stochastiques
    success = 0 
    for i in 1:100
        sim = simulation(T, s; stochastic=true, generations=200)
        final = sim[:, end] 
        vegetation = final[2] + final[3] + final[4] # quantité de végétation totale
        shrubs = final[3] + final[4] # quantité de buissons totale
        
        # conditions pour que la simulation stochastique soit un succès
        cond1 = abs(vegetation-40) <= 5 # 20% de végétation 
        cond2 = abs(final[2] - 12) <= 5 # 6% d'herbe 
        cond3 = min(final[3], final[4]) >= 0.3*shrubs # 30% du buisson le moins abondant
        if cond1 && cond2 && cond3 
        success += 1
        end
     end
     return success/100 
end

println("Success rate = ", check_success(T, s)) # on montre la fraction de succès des simulations stochastiques
# # Discussion
# Nous avons appliqué la simulation à un corridor de 200 parcelles, dont une partie des parcelles est initialement plantée avec deux espèces de buissons pour 
# simuler un aménagement du corridor sous une ligne électrique à haute tension. Nos résultats montrent que le système atteint un équilibre relativement 
# rapidement, après environ 15 générations. Nous avons réussi à simuler une gestion du corridor qui vise à limiter la croissance de végétation dense sous la ligne 
# électrique en obtenant une majorité des parcelles vides à l'équilibre. Une proportion plus faible de parcelles est occupée par de la végétation, l’herbe et les 
# buissons. Les proportions observées à l’équilibre respectent globalement les objectifs du mandat, soit 19,92 % de parcelles végétalisées, dont 29,52 % d’herbes 
# et 70,48 % de buissons. En revanche, nous n'avons pas obtenu les bonnes proportions relatives des buissons. Celles-ci étaient relativement proches, le buisson 
# le moins important représentant 48,80 % des parcelles occupées par les buissons. 
# Les simulations stochastiques montrent une variabilité autour de la trajectoire déterministe car elles reflètent le caractère aléatoire des transitions entre 
# états dans un système écologique réel. Malgré cela, les proportions à l'équilibre restent relativement stables, indiquant que la stratégie de plantation et la 
# matrice de transition choisies permettent d’atteindre un équilibre satisfaisant entre biodiversité et gestion du corridor. Ce modèle reste néanmoins simpliste 
# car il ne prend pas en compte d’autres facteurs écologiques l'influence de l'abondance de la végétation sur la croissance des plantes ou encore les perturbations 
# environnementales.


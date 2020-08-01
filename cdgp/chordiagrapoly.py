# ############ Curves, Knots,, Chord Diagrams, Graphs, Polynomials ############

"""
This subsubmodule contains functions dealing with : 
    Chord diagrams of loops in the plane (for instance coming from knot diagrams)
    Interlace graphs of chord diagrams (with orientations and signs)
    Graph labeled tree factorisation (aka Cunningham's modular decomposition)
    Fast computations for energy partition functions

    The aim is to understand the relationships betwen various polynommial invariants 
    associated to those topological and combinatorial objects such as
    Fricke polynomials of modular geodesics, Kauffman brackets of modular knots
    Energy partition functions of chord diagrams, Tutte polynomials of interlace graphs

    Main applications in mind:
    Computing partition functions of chord diagramms comming from knotted DNA
    Arithmetics of Gauss class groups and sandpile groups
"""


### ### LIBRARIES ### ###

from typing import List, Callable  
import numpy as np
import cypari
#import matplotlib.pyplot as plt


## ###### Classes, choix des structures de données, et conversions ######

class ChorDiag():
    """
    :code: ChorDiag(['a', 'b', 'c', 'd', 'b', 'a', 'd', 'c'],
                    orientations = [1,0,1,0,1,0,1,0],
                    signs = [-1,1,-1,1,-1,1,-1,1])
    *Args:
        'labels' (list) : list of even length, the labels, each label appearing twice
        'orientations' (list of bool) : list of orientations first vector second vector
        'signs'` (list of bool) : list of signs 
    """

class Graph:
    

## ###### FROM LYNDON WORDS TO CHORD DIAGRAMS ######

"""
FONCTION 2 : De la courbe au diagramme de cordes

ENTREE : Mot L&R)
SORTIE : son diagramme de cordes (linéaire), c'est à dire un mot dans lequel chaque symbole apparait exactement deux fois par exemple abcbdadc
PROCEDE : 
on parcoure la courbe paramétréé, 
à chaque fois qu'on rencontre une intersection :
    soit elle n'a pas de label, on lui en donne un et on l'append au mot, on continue
    soit elle a déjà déjà un label, on l'append au mot, on continue
arrivée au point de départ on a notre diagramme de cordes
"""

def chordiag_of_lyndon(word : str) -> list of list:


## ###### FROM CHORD DIAGRAMS TO INTERLACE GRAPHS ######

"""
FONCTION 3 : Du diagramme de cordes au graphe d'enlacement

ENTREE : un diagramme de cordes
SORTIE : un graphe
PROCEDE : 
pour chaque paire de lettres semblable (same-label) il y a un sommet, 
deux lettre sont reliées par une arête lorsqu'elles s'enlacent cycliquemenet
par exemple ..a...b...b...a.. ne s'enlacent pas mais ..a...b...a...b.. s'enlacent
"""

def gradma_of_chordiag(chordiag : list of str) -> np.matrix of bool :

    """
    Returns the interlace graph of an oriented (signed or not) chord diagram
    in the form of its adjacency matrix (gradma for graph adjaency matrix)
    The entries are :
    -1 for two  0, 1, 

    :code:`gradma_of_chordiag(ChorDiag([['a', 1], ['b', 0], ['c', 1],  ['d', 0],
                                        ['b', 1], ['a', 0], ['d', 1], ['c', 0]])

    Args:
        `chordiag` (ChorDiag): an oriented chord diagram, no signs needed

    Returns:
        Graph: an instance of our Graph class with the adjacency matrix of the chord diagram

    :Example:
        >>> gradma_of_chordiag(Chordiag())
        'Graph()'
    """

def chordiag_genus(chordiag : list of str) -> int:
    adj_mod_2 = gradma_of_chordiag(chordiag)%2 # adjacency matrix mod 2
    rank = np.rank(adj_mod_2)
    return rank

def chordiag_is_gauss(chordiag : list of str) -> bool:




## ###### CUNNINGHAM GRAPH LABELED TREE FACTORISATION OF GRAPHS ######

"""
FONCTION 4 : Décomposition de Cunningham du graphes

ENTREE : un graphe
SORTIE : un arbre de graphes "premiers"
PROCEDE : c'est un algo pénible avec des "split", 
il faudrait se servir d'un truc déjà implémenté

"""

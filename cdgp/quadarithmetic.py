
# ############ QUADRATIC ARITHMETIC : COMPOSITION AND LINKING ############


"""
    This submodule contains functions dealing with :
        Arithmetic of quadratic numbers, forms and fields
        Linking of modular knots and intersections of modular geodesics
        Fricke polynomial
"""


### ### LIBRARIES ### ###

from typing import List, Callable, Tuple 
from math import gcd
from functools import reduce
import numpy as np
import cypari
import matplotlib.pyplot as plt


## ############ DICTIONARY : LYNDON WORDS, MATRICES AND QUADRATIC FORMS ############ ##

"""
    PARI GP deals efficiently with quadratic forms, fields, Gauss composition, et.
    Quadratic forms which are indefinite correspond to hyperbolic elements in PSL_2(Z)
    Their classes (change variable, conjugacy) are represented by cyclic binary words
    Lyndon words choose the minimal representative among circular shifts.
    Thus Lyndon words uniqueley corresponds to classes of indefinite quadratic forms.

    A couple of functions give translation between Lyndon words and quadratic forms.
    (The intermediary step is through matrices in sl2 and SL2, 
    we need continued fractions to solve a Pell equation)
    We recast some computations from PARI GP in terms of Lyndon words.
    The aim is to provide a florilège of class groups with positive discriminant.
"""


def list_of_circular_shifts(word : str) -> List[str]:

    """
    Returns the list of circular shifts of a string :
    :code:`list_of_circular_shifts('RLR')`.

    Args:
        `word` (str): a word (supposedly binary L&R)

    Returns:
        list of str: list of circular shifts

    :Example:
        >>> list_of_circular_shifts('RLR')
        ['RLR', 'LRR', 'RRL']
        >>> list_of_circular_shifts('LL')
        ['LL']
    """

    n = len(word)
    return ["".join([word[i - j] for i in range(n)]) for j in range(n)]


def lyndon_of_word(word : str, comp: Callable[[List[str]],str] = min ) -> str:
    """
    Returns the Lyndon representative among set of circular shifts, 
    that is the minimum for th lexicographic order 'L'<'R'

    :code:`lyndon_of_word('RLR')`.

    Args:
        `word` (str): a word (supposedly binary L&R)
        `comp` ( Callable[List[str],str] ): comparision function  min or max

    Returns:
        str: list of circular shifts

    :Example:
        >>> lyndon_of_word('LRRLRLL')
        'LLLRRLR'
   """

    if word == '':
        return ''
    return comp(list_of_circular_shifts(word))


def compress(word : str ) -> str :

    """
    Returns compact representation of word.
    :code:`compress('LLRRRR')`

    Args:
        `word` (str): a word (supposedly binary L&R)

    Returns:
        list of str: list of circular shifts

    :Example: 
        >>> compres('LLRRRR')
        L2 R4
    """
 
    letters = ['L', 'R']
    first_letter = word[0]
    current_letter = word[0]
    counts = [0]
    for c in word:
        if c == current_letter:
            counts[-1] += 1
        else:
            current_letter = c
            counts.append(1)
 
    compress = ""
    for i,c in enumerate(counts):
        choose = i%2 if first_letter == 'L' else 1-i%2
        compress += letters[choose]+ ("" if c == 1 else str(c)) + " "
 
    return compress[:-1]


def matrix_of_word(word : str) -> np.array:
    
    """
    Returns the SL_2(Z) matrix corresponding to a given 'L'/'R' address.
    :code:`matrix_of_word('LLRR')`

    Args:
        `word` (str): a word (supposedly binary L&R)

    Returns:
        `prod` (np.array) : a numpy 2 by 2 matrix

    :Example: 
        >>> matrix_of_word('LLRR')
        array([[1, 2], [2, 5]])
    """

    R = np.array([[1,1],[0,1]])
    L = np.array([[1,0],[1,1]])
    prod = np.array([[1,0],[0,1]])
    for c in word:
        if c == 'R':
            prod = prod.dot(R)
        elif c=='L':
            prod = prod.dot(L)
    return prod
 

def bqf_of_word(word : str) -> Tuple[int]:
    
    """
    Returns the triple l,m,r) representing the quadratic form 'lx^2+m*x*y*y^2'
    corresponding to a given 'L'/'R' address.

    :code:`bfq_of_word('LLRR')`

    Args:
        `word` (str): a word (supposedly binary L&R)

    Returns:
        `(l,m,r)` (tuple of int) : three integer coefficients

    :Example: 
        >>> bfq_of_word('LLRR')
        (1,1,-1)
    """

    A = matrix_of_word(word)
    a,b,c = A[1,0], A[1,1]-A[0,0], -A[0,1]
    u = reduce(gcd,[a,b,c])
    l,m,r = a//u, b//u, c//u
    return(l,m,r)


# A paritr de là c'est du chantier

def word_of_bfq(l,m,r):
    """
    Returns the SL_2(Z) matrix corresponding to a given 'L'/'R' address.
    It works like this :
    - firsts choose the Lyndon rpresentative of word
    - computes its trace to get discriminant D of quadratic field (check it is fundamental) 
    - set w = ((0 or 1)+sqrt(D))/2 and compute fundamental unit epsilon = t + u sqrt(D)
    - and the associated fundamntal unit of the corresponding quadratic number field.

    :code:`matrix_of_word('LLRRRR')`

    Args:
        `word` (str): a word (supposedly binary L&R)

    Returns:


    :Example: 
        >>> compres('LLRRRR')
        L2 R4
    """
    
    epsilon = pari.quadunit(D)
    t = pari.trace(epsilon)
    u1 = (epsilon-t)/2
    u2 = (u1-pari.conj(u1))/2
    u = pari.norm(v2)/D

    return word


### From Lyndon words to modular geodesics ###

"""
Cette section ## n'est pas une priorité tant que je ne suis pas certain qu'elle soit incontournable. 
Ca peut être amusant et agréable de dessiner la courbe "tortue" à partir du mot de Lyondon
mais ce ne sera pas facile d'en déduire le diagramme de cordes, 
et je pense qu'on peut trouver une manière de courcirctuier le dessin 
et déduire le diagramme de cordes à partir du mot de Lyndon : il faut que je réfléchisse ! 
"""

def loop_of_lyndon(word : str) -> None:
    """
    FONCTION 1 : Du mot à la Courbe
    ENTREE : un mot en L et R, du genre LLLRRLRLLRLLRR
    SORTIE : la praramétrization d'un chemin dans le plan privé de deux points
    (je pense à un chemin du type "tortue" comme pour la fourmi de Langton, 
    ou bien à une suite de segments avec des virages à gauche ou à droite)
    PROCEDE : il faudrait que je te détaille à l'oral comment construire le chemin
    """
    pass

## ######################## RADMACHER AND LINKING ######################## ##

def Rademacher(word):
    return word.count("L")-word.count("T")

### L'algo de Pierre pour le crossing number ###

def orderlex(l1, l2, n1, n2):
    """ lexicographic order on infinite words, returns : 
    0 if the n1-th shift of l1^infty is smaller than the n2-th shift of l2^infty, 1 if it is larger, and 2 if the two words coincide 
    (this can be determined by looking only length(l1)+length(l2) letters since u^infty = v^infty iff uv = vu).
    """
    i = 0
    while( (l1[(i+n1)%len(l1)] == l2[(i+n2)%len(l2)]) & (i <= (len(l1)+len(l2))) ):
        i = i+1
    if l1[(i+n1)%len(l1)] < l2[(i+n2)%len(l2)]:
        return 0
    else:
        if l1[(i+n1)%len(l1)] > l2[(i+n2)%len(l2)]:
            return 1
        else:
            return 2

def cross_word(w1, w2):
    """ The function cross computes the crossing number of l1 and l2 on the template.
    It relies on the observation that a crossing occurs 
    when an arc coming from the left ear (the 0) goes to the right 
    of an arc coming from the right ear (the 1).
    """
    c = 0
    for i in range(len(w1)):
        for j in range(len(w2)):
            if ( (w1[i] == 'L') & (w2[j] == 'R') & (orderlex(w1,w2,i+1,j+1) == 1) ):
                c = c+1
            if ( (w1[i] == 'R') & (w2[j] == 'L') & (orderlex(w1,w2,i+1,j+1) == 0) ):
                c = c+1
    return c

### Ma formule combinatoire, algorithmiquement ###

""" We check that cross_word = enl on systematic samples of examples :
    it is always true : youpii !!
    cross_word is much faster so must be usd algoruthmically
    although enl has mathematical advantages as it consists in a closd formula

"""


def linking_patterns(max_word_length, current_pair=('LR', 'RL')):
    """ Returns the word base for the computation of enl.
    
    # The set of pairs form a binary tree which we fill in à la Pascal
    
    # The root pair is ('LR','TR')
    # The pairs to the extreme left are ('Ln R', 'R Ln')
    # The pairs on the extreme right are ('L Rn', 'Rn L')

    # The children of (G, D) are
    # to the left (LP, LQ) with LP=G[:-1]+'LR' (on enlève 'R' on rajoute 'LR') and LQ=D+'L'
    # to the right (RP, RQ) with RP=G+'R' and RQ=D[:-1]+'RL' (on enlève 'L' on rajoute 'RL')
    
    # Note the properties : 'G' ends with 'T' and 'D' ends with 'L' are preserved
    # which is why G[:-1] = G - 'T' and D[:-1] = D - 'L'
    """
 
    if len(current_pair[0]) > max_word_length:
        return []
 
    pair_left = current_pair[0][:-1]+'LR', current_pair[1]+'L'
    pair_right = current_pair[0]+'T', current_pair[1][:-1]+'RL'
 
    return [current_pair] + linking_patterns(max_word_length, pair_left) +\
                            linking_patterns(max_word_length, pair_right)
 
def occ(P,A):
    # Returns the number of times P appears at the begining of circular shifts of A.
    shifts = list_of_circular_shifts(A)
    counter = 0
    n,l = len(P), len(A)
 
    for shift in shifts:
        power = shift + shift * (n//l)
        if power[:n] == P:
            counter += 1
 
    return counter
 
def scal_PQ(P,Q,A,B):
    return (occ(P,A)*occ(Q,B)+ occ(P,B)*occ(Q,A))
 
def enl(A,B):
    # Returns the enl metric on the words A and B in the L/R alphabet.
    patterns = linking_patterns(len(A)+len(B)+1)
    return sum([scal_PQ(P,Q,A,B) for P,Q in patterns])


## ###### FIRCKE POLYNOMIAL (TRACE OF QUANTUM DEFORMATION) ######


def quantize(word):
    Lq = cypari.pari("[q,0;1,1/q]")
    Tq = cypari.pari("[q,1;0,1/q]")

    if len(word) == 0:
        return cypari.pari("[1,0;0,1]")
    result = Lq if word[0] == "L" else Tq

    for c in word[1:]:
        result *= (Lq if c == "L" else Tq)

    return result 

def Fricke(word):
    M = quantize(word)
    return M[0][0]+M[1][1]

def are_Fricke_equiv(word1,word2):
    return Fricke(word1) == Fricke(word2)    


## ###### FLORILEGE OF CLASS GROUPS ######




################################################ PIECES EN STOCK ################################################


### Discriminant, resultant and Killing form

def discriminant(mat):
    return (mat.trace()**2-4)

def killing_form(matrix_A,matrix_B):
    #m_A = matrix_of_address(word_a) 
    #m_B = matrix_of_address(word_b)
    scal = 2* np.matmul(matrix_A, matrix_B).trace()-matrix_A.trace()*matrix_B.trace()
    return scal

def resultant(mat_A,mat_B):
    la,ma,ua = quadratic_form_of_matrix(mat_A)
    lb,mb,ub = quadratic_form_of_matrix(mat_B)
    res = (la*ua-lb*ub)**2-ma*mb*(la*ua+lb*ub)+(la*mb**2*ua)+(lb*ma**2*ub)
    return res


## INVERSE DE GAUSS

def Gauss_inverse(word):
    return word[::-1]

def is_ambiguous(word):
    return circular_min_rpz(word) == circular_min_rpz(Gauss_inverse(word))

def one_complent(word):
    to_return = ""
    for c in word:
        if c == 'L':
            to_return += 'T'
        else:
            to_return += 'L'

    return to_return

def is_inert(word):
    return circular_min_rpz(word) == circular_min_rpz(one_complent(word))

def is_reciprocal(word):
    return circular_min_rpz(word) == circular_min_rpz(Gauss_inverse(one_complent(word)))

## GENUS

def is_principal_genus(word):
    rep = circular_min_rpz(word)
    return (compress(rep).count('L') == 1) and (compress(rep).count('T') == 1)

## INVERSE MATRICIEL

def matrix_inverse(mat):
    a,b = mat[0][0],mat[0][1]
    c,d = mat[1][0],mat[1][1]
    return np.array([[d,-b],[-c,a]])



## FRACTIONS CONTIUES


def continued_fraction(x,long, expansion=[], A = np.matrix([[1,0],[0,1]]), colonnes=[] ):
    """ 
    entrer : un réel x et un entier long 
    retourne :
    expansion : 
    colonnes : liste des 'long' premières approximations par fractions continuées partielles de 'x'
    """
    
    if long <= 0:
        return(expansion, colonnes, A)
    
    n=int(np.floor(x))
    expansion.append(n)
    
    A=np.matmul(A,np.matrix([[n,1],[1,0]]))
    colonnes.append("{} /{}".format(int(A[0,1]), int(A[1,1]) ))
    
    return continued_fraction(1/(x-n),long-1, expansion, A, colonnes)


def continued_fraction_hj(x,long, expansion=[], A = np.matrix([[1,0],[0,1]]), colonnes=[] ):
    """ Idem mais pour les fractions continues d'Hirzebruch Youg
    """
    if long <= 0:
        return expansion#, colonnes, A
    
    n=int(np.ceil(x))
    expansion.append(n)
    
    A=np.matmul(A,np.matrix([[n,-1],[1,0]]))
    colonnes.append("{} /{}".format(int(A[0,1]), int(A[1,1]) ))
    
    return continued_fraction_hj(1/(n-x),long-1, expansion, A, colonnes)



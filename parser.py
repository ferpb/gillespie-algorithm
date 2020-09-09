# -----------------------------------------------------------------------------
# calc.py
#
# A simple calculator with variables.
# -----------------------------------------------------------------------------

################################
#            LEXER             #
################################

tokens = [
    'INTEGER',
    'FLOAT',
    'NAME',
    'PLUS',
    'ARROW',
    'EQUAL',
    'SEMICOLON',
    'COLON'
]

# Tokens

t_PLUS      = r'\+'
t_EQUAL     = r'='
t_ARROW     = r'->'
t_SEMICOLON = r';'
t_COLON     = r':'

# Reserved words
reserved = {
    'initial_time' : 'INITIALTIME',
    'final_time' : 'FINALTIME',
}

tokens += list(reserved.values())

def t_NAME(t):
    r'[a-zA-Z_][a-zA-Z0-9_]*'
    t.type = reserved.get(t.value,'NAME') # check for reserved words
    return t

def t_FLOAT(t):
    r'\d+\.\d+'
    t.value = float(t.value)
    return t

def t_INTEGER(t):
    r'\d+'
    t.value = int(t.value)
    return t


# Ignored characters
t_ignore = " \t"
t_ignore_COMMENT = r'\#.*'

def t_newline(t):
    r'\n+'
    t.lexer.lineno += t.value.count("\n")

def t_error(t):
    print(f"Illegal character {t.value[0]!r}")
    t.lexer.skip(1)


# Build the lexer
import ply.lex as lex
lexer = lex.lex()


# lexer.input(data)

# while True:
#     tok = lexer.token()
#     if not tok:
#         break
#     print(tok)


################################
#           PARSER             #
################################

# Dictionary of reactions
# reaction format: {<name>: {'reactives': {<name>: <quantity>,...},
#                  {'products': {<name>: <quantity>,...},
#                  {'constant': <name>}}
reactions = {}
# Dictionary of names for storings species concentrations and constants values
# name format: {<name>: <value>}
names = {}
initial_time = 0
final_time = 100

def p_data(p):
    'data : initialtime finaltime reaction_list assignation_list'

def p_reaction_list(p):
    '''reaction_list : reaction_list reaction
                     | reaction'''

def p_reaction(p):
    'reaction : NAME COLON species_list ARROW species_list SEMICOLON NAME'
    global reactions
    reactions[p[1]] = {'reactives': p[3], 'products': p[5], 'constant': p[7]}

def p_species_list(p):
    '''species_list : species_list PLUS specie
                    | specie '''
    species = {}
    if (len(p) > 2):
        species = p[1]
        (specie, num) = p[3]
        species[specie] = num
    else:
        (specie, num) = p[1]
        species[specie] = num

    p[0] = species

# return (name, num)
def p_specie(p):
    '''specie : NAME
              | INTEGER NAME'''
    if (len(p) > 2):
        p[0] = (p[2], p[1])
    else:
        p[0] = (p[1], 1)

def p_quantity(p):
    '''quantity : INTEGER
                | FLOAT'''
    p[0] = p[1]

def p_initialtime(p):
    'initialtime : INITIALTIME EQUAL quantity'
    global initial_time
    initial_time = p[3]

def p_finaltime(p):
    'finaltime : FINALTIME EQUAL quantity'
    global final_time
    final_time = p[3]

def p_assignation_list(p):
   '''assignation_list : assignation_list assignation
                       | assignation'''

def p_assignation(p):
    'assignation : NAME EQUAL quantity'
    global names
    names[p[1]] = p[3]

def p_error(p):
    print(f"Syntax error at {p.value!r}")


import ply.yacc as yacc
yacc.yacc()

def parse(data):
    yacc.parse(data)
    return (initial_time, final_time, reactions, names)



################################
#      GILLESPIE PARSER        #
################################

import numpy as np

# NOTE: From Python 3.7 the built in dictionary is guaranteed to maintain
# the insertion order

def gillespie_parse(data):
    (initial_time, final_time, reactions, names) = parse(data)
    print(initial_time)
    print(final_time)
    print(reactions)
    print(names)
    print()
    print("Empezar gillespie")

    # Get all the species in the instance
    species = []
    for rkey, rvalue in reactions.items():
        print(rkey, rvalue)
        for skey in rvalue['reactives']:
            print(skey)
            species.append(skey)
        for skey in rvalue['products']:
            print(skey)
            species.append(skey)
    # Remove duplicates from species list maintaining order
    species = list(dict.fromkeys(species))
    print(species)

    # Construct update matrix (num_reactions x num_species)
    # and propensity_functions list
    update_matrix = np.zeros((len(reactions), len(species)))
    propensity_functions = []

    for i, rvalue in enumerate(reactions.values()):
        print(rvalue)
        for skey, svalue in rvalue['reactives'].items():
            # reactive (negative)
            update_matrix[i][species.index(skey)] = -svalue
        for skey, svalue in rvalue['products'].items():
            # product (positive)
            update_matrix[i][species.index(skey)] = svalue

        try:
            propensity_functions.append(names[rvalue['constant']])
        except KeyError:
            propensity_functions.append(0)

    # Construct species concentration list
    species_concentration = []
    for s in species:
        try:
            species_concentration.append(names[s])
        except KeyError:
            species_concentration.append(0)

    print(initial_time)
    print(final_time)
    print(update_matrix)
    print(species)
    print(propensity_functions)
    print(species_concentration)
    return(initial_time, final_time, update_matrix, species, propensity_functions, species_concentration)

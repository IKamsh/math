def is_solution_finite(system):
    '''
    'system' is expected to be a list of polynomials
    
    the function does not take into account the case
    when system has no solution
    '''
    
    R = sum(system).parent()
    variables = R.gens() # variables array
    n = len(variables) # number of variables
    I = ideal(system)
    J = I.groebner_basis()
    
    if n > len(J):
        return False
    
    for i in range(n): # var cycle
        for j in range(len(J)): # polynom cycle
            counter = 0
            if J[j].lm() % variables[i] == 0:
                counter += 1
                for k in range(n): # check that other variables do not divide
                    if k != i:
                        if J[j].lm() % variables[k] == 0:
                            counter = 0
                            break
            
            if counter == 1: #found proper polynom and move to the next variable
                break
            if j == len(J) - 1:  # if no proper polynoms
                return False
                            
    return True   
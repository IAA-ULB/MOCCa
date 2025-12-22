import itertools

def enumerate_index_combinations(ndim, arg, coupling):
  """
  Enumerate all possible combinations of indices for a tensor contraction. 

  Examples:

    1.If you try to get the indices for the following calculation of a density:
  
        D^N(N,N)N_mn =  sum_k D^NN,NN_mkkn 

      You can call this routine with 
         - ndim = 4, i.e. the TOTAL number of indices of the lhs tensor
         - args = (m,n)
         - coupling = [(1,2)], i.e. a coupling between the second and third index of the lhs tensor

      And it should result in 

        sum_args = [(m,0,0,n),(m,1,1,n),(m,2,2,n)]
        
      i.e. the relevant indices to explicitly code for the sum in the rhs.

  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Input: 
    ndim   : int
      total number of dimensions of the lhs tensor (i.e. including contracted indices).
    arg    : tuple
      indices of the left-hand side tensor, i.e. concrete values for the uncontracted indices. 
      assumed to be in the ordering of the indices of the lhs tensor.

    coupling: list of tuples
      a list of two-index couplings, i.e. which indices of the lhs tensor are coupled/contracted.
      
  Output:
    sum_args: list of tuples
      a list of all possible combinations of indices for the rhs tensor, given the couplings.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  """

  Nc = len(coupling)

  sum_args = []
  if Nc == 0:
    # Nothing to do if no couplings needed
    uncontracted = [arg]
  else:
    # These are all possible combinations for the couplings
    cont = itertools.product(range(3), repeat=Nc)

    # For each possible value of the couplings, build the full argument tuple full_arg
    for c in cont:                          
      ii = 0
      full_arg = ()

      # For each index, check if it is coupled or not
      for index in range(ndim):              
        coupled = False 
        for ic, coup in enumerate(coupling):
          if(index in coup):
            coupled = True
            break
        
        if(not coupled):
          # Not coupled -> use the value of a free index from arg
          full_arg += (arg[ii],)
          ii += 1
        else:
          # Coupled -> put the corresponding value from c
          full_arg += (c[ic],)
      sum_args.append(full_arg)
  return sum_args

def reduce_vector_coupling(ndim, coupling, cross):
  """
    Transform a complicated problem featuring 
      - 2-index couplings
      - 3-index couplings (cross-couplings that are summed over)
      - vector products (cross-couplings with a free index)

    into a simpler problem featuring only 2-index couplings
    by adding explicit Levi-Civita indices to problem. 

    Example:

      An expression where a vector product is summed over 
          Derm_D_I_I D_I_NxmSxm 
      should have as input  
          ndim     = 3 
          coupling = [(0,1,2)]
          cross    = []
      and should result in 
          new_coupling = [(0,3),(1,4),(2,5)]
      which corresponds to an object that effectively looks like this
          Derm_D_I_I D_I_NiSj E_mij 
      where E_mij is the Levi-Civita symbol.

    - - - - - - - - - - - - - - - - - - - - - - -
    Input:
    ndim   : int
      total number of indices of the uncontracted tensor 
    coupling: list of tuples
      a list of couplings, each being either a 2-index coupling 
      or a 3-index coupling (vector product summed over)
    cross   : list of tuples
      a list of cross-couplings, each being a 2-index coupling 
      representing a vector product with a free index.
    - - - - - - - - - - - - - - - - - - - - - - -
    Output: 
    new_ndim   : int
      total number of indices of the uncontracted tensor after 
      possible addition of Levi-Civita indices.
    new_coupling: list of tuples
      a list of couplings, each being a 2-index coupling only.
    N_levicivita: int
      number of Levi-Civita indices added to the problem.
    - - - - - - - - - - - - - - - - - - - - - - -
  """

  new_ndim = ndim
  new_coupling = []
  N_levicivita = 0
  for c in coupling:
    if(len(c) == 2):
      # 2-length couplings need no further work
      new_coupling.append(c)
    elif(len(c) == 3):
      # 3-length couplings need to be expanded
      new_coupling.append( (c[0], new_ndim) )
      new_coupling.append( (c[1], new_ndim+1) )
      new_coupling.append( (c[2], new_ndim+2) )
      new_ndim += 3
      N_levicivita += 1

  for c in cross:
    new_coupling.append( (c[0], new_ndim+1) )
    new_coupling.append( (c[1], new_ndim+2) )
    new_ndim += 3
    N_levicivita += 1

  return new_ndim, new_coupling, N_levicivita

def recombine_levicivita(args, N_levicivita):
  """
    Convert a list of indices for object with virtual Levi-Civita symbols added 
      back into a list of indices for objects without Levi-Civita symbols with 
      additional signs.

    See contract_levi_civita_indices for the operation done on each individual tuple.

  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Input: 
    args : list of tuples
      The list of argument tuples to be processed.
    N_levicivita : int
      The number of Levi-Civita symbols to be contracted.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Output:    
    actual_args : list of tuples
      The list of argument tuples after contracting the Levi-Civita indices.
    signs : list of int
      The list of signs resulting from the contraction of the Levi-Civita symbols.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  """
  actual_args = []
  signs       = []
  for arg in args:
    actual_arg, sign = contract_levi_civita_indices(arg, N_levicivita)
    
    if(sign != 0):
      actual_args.append(actual_arg)
      signs.append(sign)

  return actual_args, signs

def contract_levi_civita_indices(arg, N_levicivita):
  """
    Turn a tuple of indices for an object with virtual Levi-Civita symbols added
    into a tuple of indices for an object without Levi-Civita symbols with
    an additional sign.

    Example: 
      An expression where a vector product is summed over with a Levi-Civita symbol

          Der_D_I_NS E with argument tuple (0,1,2,0,1,2)
                           | 
                           -> 1 added/virtual Levi-Civita symbol

      should result in the arguments
          (0,1,2) with sign +1
          (0,2,1) with sign -1
      that label the abstract object 
          Der_D_I_NS
      and - when summed - result in 
          Derm_D_I_NxmSxm

  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Input: 
    arg: tuple
      The argument tuple containing indices including Levi-Civita indices.
    N_levicivita: int
      The number of virtual Levi-Civita symbols that were added.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Output:
    actual_arg: tuple
      The argument tuple after contracting the Levi-Civita indices.
    sign: int
      The sign resulting from the contraction of the Levi-Civita symbols.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
  """

  sign = 1
  if(N_levicivita == 0):
    return arg, sign
  
  # The last (3 * N_levicivita) indices are the Levi-Civita indices
  levi_list = list(arg[-3 * N_levicivita:]) 
  if len(set(levi_list)) < len(levi_list):
    # If there are repeated indices, the Levi-Civita symbol is zero
    sign = 0
    actual_arg = []
  else:
    # The first len(arg)-3 * N_levicivita indices are the actual indices
    actual_arg = arg[:-3 * N_levicivita]
    # Calculate the sign based on the permutation of the Levi-Civita indices
    for k in range(N_levicivita):
      levi_indices = levi_list[3*k:3*(k+1)]
      # Count the number of inversions to determine the sign
      for i in range(3):
        for j in range(i+1, 3):
          if levi_indices[i] > levi_indices[j]:
            sign *= -1

  return actual_arg, sign

if( __name__ == "__main__" ):

  print (" Testing a term in the U potential: 'Derm D_I_NxmSxm'")
  coupling = [(0,1,2)]
  cross    = []
  new_ndim, new_coupling, N_levicivita = reduce_vector_coupling(3, coupling, cross)
  print("Dim of new object = ", new_ndim, ' with coupling ', new_coupling, ' and N_levicivita ', N_levicivita)

  args = enumerate_index_combinations(new_ndim, (), new_coupling)
  print ("Number of terms to sum over:", len(args))
  print ("Args", args)

  actual_args, signs = recombine_levicivita(args, 1)
  print ("Number of actual terms to sum over:", len(actual_args))
  for i, actual_arg in enumerate(actual_args):
    print ("  Actual arg:", actual_arg, " with sign ", signs[i])

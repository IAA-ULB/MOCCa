import itertools

def build_index_combinations(ndim, nfree, coupling, cross):
  """
  Build all possible combinations of indices for a tensor contraction involving both standard 
  couplings and vector products.

  Examples:

    1. If you try to get the indices for the following calculation of a density:
        D^N(N,N)N_mn =  sum_k D^NN,NN_mkkn
      You can call this routine with
          - ndim = 4, i.e. the TOTAL number of indices of the lhs tensor
          - nfree = 2, i.e. the number of free indices (m,n)
          - coupling = [(1,2)], i.e. a coupling between the second and third index of the lhs tensor
      And it should result in
        {
          (0,0): [(0,0,0,0),(0,1,1,0),(0,2,2,0)],
          (0,1): [(0,0,0,1),(0,1,1,1),(0,2,2,1)],
          (0,2): [(0,0,0,2),(0,1,1,2),(0,2,2,2)],
          (1,0): [(1,0,0,0),(1,1,1,0),(1,2,2,0)],
          (1,1): [(1,0,0,1),(1,1,1,1),(1,2,2,1)],
          (1,2): [(1,0,0,2),(1,1,1,2),(1,2,2,2)],
          (2,0): [(2,0,0,0),(2,1,1,0),(2,2,2,0)],
          (2,1): [(2,0,0,1),(2,1,1,1),(2,2,2,1)],
          (2,2): [(2,0,0,2),(2,1,1,2),(2,2,2,2)]
        }
      Each entry in the dictionary corresponds to one equation that involves multiple terms. 
      For example,  Key (0,0) with value [(0,0,0,0),(0,1,1,0),(0,2,2,0)] stands for 
        D^N(N,N)N_00 = D^NN,NN_0000 + D^NN,NN_0110 + D^NN,NN_0220  
      
  Input:
    ndim   : int
      total number of dimensions of the lhs tensor (i.e. including contracted indices).
    nfree  : int
      number of free (i.e. uncontracted) indices of the lhs tensor.
    coupling: list of tuples
      a list of two-index or three-index couplings, i.e. which indices of the lhs tensor are coupled/contracted.
    cross   : list of tuples
      a list of two-index cross-couplings, i.e. which indices of the lhs tensor are coupled via a vector product.

  Output:
    dict_args: dict
      a dictionary where keys are tuples of free indices and values are lists of tuples
      representing all possible combinations of indices for the rhs tensor, given the couplings.

  """

  # All possible combinations for the free indices, i.e. the uncontracted ones
  args = itertools.product(range(3), repeat=nfree)

  # Transform the problem to a simpler one with only 2-index couplings
  new_ndim, new_coupling, N_levicivita = reduce_vector_coupling(ndim, coupling, cross)

  dict_args = {}
  for arg in args:  
    # All possible combinations for the indices of the extended problem 
    extended_args = enumerate_index_combinations(new_ndim, arg, new_coupling)
    # Recombine the Levi-Civita indices to get actual indices and signs
    actual_args = recombine_levicivita(extended_args, N_levicivita)
    dict_args[arg] = actual_args 

  return dict_args

def enumerate_index_combinations(ndim, arg, coupling):
  """
  Enumerate all possible combinations of indices for a tensor contraction with 2-index couplings.

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

  Input: 
    ndim   : int
      total number of dimensions of the lhs tensor (i.e. including contracted indices).

    arg    : tuple
      indices of the left-hand side tensor, i.e. concrete values for the uncontracted indices. 
      assumed to be in the ordering of the indices of the lhs tensor.

    coupling: list of tuples
      a list of two-index couplings, i.e. which indices of the lhs tensor are coupled/contracted.
      
  Output:
    sum_args: list of lists
      a list of all possible combinations of indices for the rhs tensor, given the couplings.

    ATTENTION: the output is a list of lists, NOT a list of tuples! This is to allow mutability
               when including signs from Levi-Civita contractions later on.
  """

  Nc = len(coupling)

  sum_args = []
  if Nc == 0:
    # Nothing to do if no couplings needed
    sum_args = [list(arg)] # Ensure that we return a list of lists!
  else:
    # These are all possible combinations for the couplings
    cont = itertools.product(range(3), repeat=Nc)

    # For each possible value of the couplings, build the full argument tuple full_arg
    for c in cont:                          
      ii = 0
      full_arg = []

      # For each index, check if it is coupled or not
      for index in range(ndim):              
        coupled = False 
        for ic, coup in enumerate(coupling):
          if(index in coup):
            coupled = True
            break
        
        if(not coupled):
          # Not coupled -> use the value of a free index from arg
          full_arg += [arg[ii]]
          ii += 1
        else:
          # Coupled -> put the corresponding value from c
          full_arg += [c[ic]]
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

    Input:

      ndim   : int
        total number of indices of the uncontracted tensor 

      coupling: list of tuples
        a list of couplings, each being either a 2-index coupling 
        or a 3-index coupling (vector product summed over)

      cross   : list of tuples
        a list of cross-couplings, each being a 2-index coupling 
        representing a vector product with a free index.

    Output: 

      new_ndim   : int
        total number of indices of the uncontracted tensor after 
        possible addition of Levi-Civita indices.

      new_coupling: list of tuples
        a list of couplings, each being a 2-index coupling only.

      N_levicivita: int
        number of Levi-Civita indices added to the problem.
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
      Signs are included by multiplying the first index by a sign.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  """
  actual_args = []
  for arg in args:
    actual_arg, sign = contract_levi_civita_indices(arg, N_levicivita)
    
    if(sign != 0):
      if (len(actual_arg) > 0):
        actual_arg[0] *= sign
      actual_args.append(tuple(actual_arg))

  return actual_args

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

   Input: 

    arg: tuple
      The argument tuple containing indices including Levi-Civita indices.

    N_levicivita: int
      The number of virtual Levi-Civita symbols that were added.

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

def regroup_indices_left( args, ldim, lcoupl):
  """
    This function regroups a dictionary containing summation indices along the conventions of
    build_index_combinations into a dictionary with different conventions that are (somewhat) optimal
    to generate contributions to the single-particle hamiltonian associated with a given mean-field 
    density or potential

            D_LeftOperator_RightOperator   or   F_LeftOperator_RightOperator

    Example:
        The set of index combinations for a term in the single-particle hamiltonian corresponding to the 
        potential F_NmNm_NkNk would be generated by build_index_combinations as

        {
          (): [(0,0,0,0),(0,0,1,1),(0,0,2,2),(1,1,0,0),(1,1,1,1),(1,1,2,2),(2,2,0,0),(2,2,1,1),(2,2,2,2),
        }

        which is not practical to generate code to calculate the action of the single-particle hamiltonian.

        Rebuilding this kind of dictionary such that terms get grouped by equal indices to the left operator
        - complete indices, including resummed ones! - makes life easier to generate efficient code - it is
        the "Left" operator which is numerically hard!

        This gives the following regrouping of the input dictionary
        {
          (0,0) : [(0,0,0,0),(0,0,1,1),(0,0,2,2)] 
          (1,1) : [(1,1,0,0),(1,1,1,1),(1,1,2,2)]
          (2,2) : [(2,2,0,0),(2,2,1,1),(2,2,2,2)]
        }

        This is still not entirely optimal, as it would be nice to replace all repeated derivative indices 
        in the left operator by calls to Laplacians. We thus eliminate all repeated left indices, both 
        in the keys and the values:

        {
          () : [(0,0),(1,1),(2,2)] 
        }

        where the remaining indices ONLY index the free indices of the leftoperator and ALL of the indices 
        of the right operator. 


    ATTENTION: this function makes the assumption that ALL of the operators on the left are derivative operators, 
               implying that any coupling that is internal to the left operator can be recombined as a Laplacian!

    Input:
      args : dictionary with tuples as keys and lists of tuples as values
        A dictionary built by build_index_combinations.

      ldim :  integer 
        the TOTAL dimension of the left operator in the corresponding mean-field density/potential

      lcoupl : list of integer tuples
        list of couplings INSIDE the left-operator; i.e. NOT ALL couplings inside the whole mean-field density/potential
      
    Output:
      regrouped_args :  dictionary with tuples as keys and lists of tuples as values
        A dictionary more optimised to the generation of terms in the single-particle hamiltonian.
        Based on args, this includes
          (i)  relabeling of keys with the values of the indices of the LeftOperator
          (ii) elimination of repeated indices for the LeftOperator
  """

  arguments = []
  for key in args.keys():
    for indices in args[key]:
      arguments.append(indices)

  regrouped_args = {}
  
  for arg in arguments:
    larg = ()
    for k in range(ldim): 
      found = False 
      for c in lcoupl: 
        if k in c:
          found = True
      if not found:
        larg += (arg[k],)

    complete_arg = larg + arg[ldim:]
    try:
      regrouped_args[larg].append(complete_arg)
    except KeyError:
      regrouped_args[larg] = [complete_arg]

  return regrouped_args 

if( __name__ == "__main__" ):

  print (" Testing a term in the U potential: 'Derm D_I_NxmSxm'")
  coupling = [(0,1,2)]
  cross    = []
  new_ndim, new_coupling, N_levicivita = reduce_vector_coupling(3, coupling, cross)
  print("Dim of new object = ", new_ndim, ' with coupling ', new_coupling, ' and N_levicivita ', N_levicivita)

  args = enumerate_index_combinations(new_ndim, (), new_coupling)
  print ("Number of terms to sum over:", len(args))
  print ("Args", args)

  actual_args = recombine_levicivita(args, 1)
  print ("Number of actual terms to sum over:", len(actual_args))
  for i, actual_arg in enumerate(actual_args):
     print ("  Actual arg:", actual_arg)


  d = build_index_combinations(3, 0, coupling, cross)
  print (d)

  print ()
  d = build_index_combinations(4, 2, [(1,2)], [])
  print (d)

  print ("Testing D_I_I")
  d = build_index_combinations(0, 0, [], [])
  print (d)
# How to write a testing script?

A MOCCa testing script has the following properties:

 - It can be invoked manually (i.e. by a user) or automatically (i.e. through github actions) by typing  
   ``` bash script.sh [options] ```  
   inside its directory. All paths must be specified with respect to that level of the directory structure.
 - It returns a single exit code indicate succes (0) or failure (1).
 - It is (nearly) self-contained meaning that it does everything in the list below.
      1. Copy all relevant files (executable(s), .param, ...) to an appropriate working directory
      2. Perform the calculations
      3. Analyse the results  
Note that this list does *__not specify compilation__*: a script is free to handle compilation or assume that the relevant executables exist.
 - If it does not handle compilation, the script must the user to specify the relevant executable through a command line option. For example: the store_derivatives test is used as follows  
 ``` bash store_derivatives.sh -p BSkG3 -e BXL```
 - It must be documented thoroughly; in particular it should have 
    1. A short explanation of what the script tests: what quantity and **to which precision**.
    2. A list of dependencies; i.e. if your script relies on something more than MOCCa and bash.
    3. A clarification of the way a script should be used, including an explanation of possible inputs.
    4. An explicit listing of the right (or at least, **expected**) output. 
    5. **A reference commit**: the hash of the commit for which you obtained the expected output.
    6. Your name (and possibly email adress) to mark you as responsible for this test.

If this list seems daunting, I recommend looking into the tests that already exist and to rely on some of the boilerplate functions in ```functions.sh``` that should render setting up bash scripts less painful.

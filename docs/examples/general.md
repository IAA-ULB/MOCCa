# Examples and automated testing

 The code comes bundled with an extensive set of scripts that serve as automated 
 tests to verify new versions of the code. 
 
 These come in three flavours:

 - integration tests: tests of all parts of the code acting in tandem.
 - regression tests: tests aimed at verifying earlier bugs have not reappeared.
 - unit tests: tests of individual parts of the code. 
 
 The integration tests also serve as inspiration to new users, detailing the operation of 
 code in what are essentially production conditions.

## Integration tests 

  - ```minimal.sh```: 
     A calculation of spherical O16 with a standard Skyrme 
     parameterisation whose energy gets compared to a known value.
  - ```symmetries.sh```:
    A showcase of different symmetry compilation options that verifies that they 
    return identical results when the symmetries are conserved self-consistently.
  - ```store_derivatives.sh```:
    Verify that calculations with store_derivatives=.true. and .false. are equivalent.

## Regression tests 

 There are not yet any regression tests that have been defined. 

## Unit tests

 There are not yet any unit tests that have been defined. 

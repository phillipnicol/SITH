## Upcoming Version TBD 

**New functionality:**

  * `simulateTumorMTBP()`, where MTBP = multi-type branching process, implements the spatial growth model without the infinite alleles hypothesis. In 
  particular, the user can define the mutations, mutation rates, and selective advantages conferred to a cell acquiring each mutation. To avoid an exponential
  number of parameters, the transitions between mutations must be described as a directed acyclic graph. 

**Changes to existing functionality:**

  * The `alleles` data frame which was returned by `simulateTumor()` has been renamed to `genotypes` to be biologically correct. 

## Version 1.0.1 July 1, 2020 

This patch contains bug fixes. In addition, the dependence on R 4.0.0 is switched to 3.6.0, so that users on the old release can use the package.

**Bug fixes**:

  * An ambiguous call to the overloaded `pow()` in "post_processing.h" causes a compilation error on Solaris. 

## Version 1.0.0 June 29, 2020

Initial submission to CRAN. 

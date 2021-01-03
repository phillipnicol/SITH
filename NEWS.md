## Version 1.1.0 January 3, 2021 

This update introduces new functionality for the tumor growth simulator. 

**New functionality:**

  * A new parameter to `simulateTumor()` has been added. This allows the user to specify a custom disease model for more control over the kinds of mutations that can occur. See vignette for more details. 
  * Coverage has been added to `bulkSample()` and `randomBulkSamples()` to better simulate real sequencing data.

**Changes to existing functionality:**

  * The `alleles` data frame which was returned by `simulateTumor()` has been renamed to `genotypes` to be biologically correct. 

**Bug fixes**:

  * In `visualizeTumor`
  
**Internal changes:**

 * The `C++` code has been completely refactored for organization. 

## Version 1.0.1 July 1, 2020 

This patch contains bug fixes. In addition, the dependence on R 4.0.0 is switched to 3.6.0, so that users on the old release can use the package.

**Bug fixes**:

  * An ambiguous call to the overloaded `pow()` in "post_processing.h" causes a compilation error on Solaris. 

## Version 1.0.0 June 29, 2020

Initial submission to CRAN. 

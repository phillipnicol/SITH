## Resubmission
This is a resubmission in which I have addressed the reviewer's comments:
  * I have added the author and doi to a relevant paper describing the methods used in the package. 
  * I have replaced \dontrun{} with \donttest{} in 'SITH-package.Rd' as the example can be run but takes longer than 5 seconds. 
  * In spatial.R I have ensured that the user's par is not changed by adding on.exit() in the relevant function. 

## Testing platforms 
  * macOS 10.15
  * win-builder (release and devel)
  * Ubuntu Xenial (travis-ci) 

0 Errors 0 Warnings 2 Notes

The first note indicates that this is a first time submission. The second note flags the word "Intra-Tumor" as 
possibly mispelled. Intra-tumor (sometimes intratumor) is a technical term that is understood by cancer researchers. 

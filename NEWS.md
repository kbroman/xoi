## R/xoi revision history

### Version 0.72, 2023-03-21

- Changelog -> NEWS.md

- Fixed "uninitialized variable" warning in C code.


### Version 0.70, 2022-01-21

- Revise xoiversion() to handle a case like "0.70".


### Version 0.69-2, 2020-02-26

- Fix use of class(), using for example inherits(cross, "cross") to
  determine if an object has class "cross". Also, added an internal
  function crosstype() for grabbing the cross type.


### Version 0.68, 2019-05-21

- Small bug fix re || vs | in function for gamma probabilities

- Use markdown in documentation


### Version 0.67, 2018-03-08

- Small bug fix in simStahl.c re over-running max_nxo.

- Small changes to avoid Note about "R_registerRoutines"

- Remove calls to missing() and use NULL argument defaults instead.


### Version 0.66, 2015-10-18

- Add example in the help files for est.recrate() and recrate2scanone().


### Version 0.65, 2014-09-23

-  Revise simStahl to: (a) use different code in the case that the
   interference parameter, nu, is an integer, and (b) include the
   option obligate_chiasma=TRUE/FALSE for that integer nu case.


### Version 0.64, 2014-05-09

-  Convert to using Roxygen2 for managing documentation


### Version 0.63, 2013-11-01

-  Revised est.recrate to handle multiple chromosomes at once.


### Version 0.62, 2013-09-13

-  Fixed a big in simStahl


### Version 0.61, 2013-01-14

-  Incorporated C code for Brent_fmin from R


### Version 0.60, 2012-11-21

-  Added some error messages for fitStahl and stahlLoglik, if input is
   not a list.

-  Modified find.breaks to handle an F2 intercross.

-  Revised stahlLoglik to handle the case of an intercross.

-  stahlLoglik returns just the vector of log likelihoods (rather than
   a data.frame with the parameters as well).  The parameter values
   are saved as attributes of the result.


### Version 0.59-2, 2012-10-15

-  Removed use of R_zeroin() by incorporating code within the package
   as the function Rxoi_zeroin().

-  Reduced computation time for examples in several help files.


### Version 0.58-3, 2012-03-02

-  Revised xoiversion() to use packageVersion().


### Version 0.58-2, 2011-11-07

-  Added NAMESPACE.

-  Compressed bssbsb data.


### Version 0.57-5, 2011-07-28

-  Added functions intensity() and coincidence(), plus internal
   functions get_n_xo and identify_xo, from Il youp Kwak


### Version 0.56-4, 2011-05-09

-  Removed a few extraneous .'s from titles in help files


### Version 0.56-3, 2010-05-04

-  Removed a few unused variables from C code.


### Version 0.56-2, 2010-04-11

-  Revised chiasma() slightly, so that "Done!" is printed only if
   verbose=TRUE.


### Version 0.56-1, 2009-10-26

-  Fixed a bug in find.breaks (changing call to locate.xo to
   locateXO).

-  Fixed links in help files.


### Version 0.55-2, 2009-06-02

-  Added functions stahlLoglik and fitStahl, for calculating the log
   likelihood under the Stahl model and for obtaining MLEs and the log
   LR for testing p=0.

-  Revised convertxoloc() so that the output is a data frame.


### Version 0.54-1, 2009-05-18

-  In est.recrate, modified the default positions at which the
   recombination rate is calculated (that is, if the "pos" argument is
   not specified) so that calculations will be done at the marker
   positions and at a grid with 4 positions per Mbp.

-  Added a function recrate2scanone for converting est.recrate results
   for multiple chromosomes to the form of the output of scanone (from
   R/qtl), but I'll leave it undocumented for now.


### Version 0.53-4, 2009-04-28

-  Changed the license to GPL-3.


### Version 0.53-3, 2009-04-28

-  Commented out an empty example in a help file.


### Version 0.53-2, 2008-10-16

-  Fixed a bug in fitGamma.


### Version 0.53-1, 2008-08-25

-  Added a function est.recrate for obtaining a smoothed estimate of
   the recombination rate along a chromosome, using the cM and Mbp
   position of markers.


### Version 0.52-8, 2008-08-21

-  Fixed a slight error in the help file for fitGamma.


### Version 0.52-7, 2007-10-09

-  Fixed some slight problems in the help files.


### Version 0.52-6, 2007-09-20

-   Changed my email address, etc., due to my move to UW-Madison


### Version 0.52-5, 2007-04-26

-   Added a function joint.given.two for calculating the joint density
    of the crossover locations, given that there are two.

-   Added a function est.coi, for estimating the three-point
    coincidence along a chromosome.


### Version 0.51-3, 2006-12-20

-   Added function est.intensity and est.kfunc to estimate the
    intensity and K_inhom functions for an inhomogeneous point
    process.

-   Revised most of the functions for calculating densities and
    such so that you may specify (via the argument "x") the
    exact values at which these functions are to be calculated.
    If x is unspecified, we calculate things at n equally spaced
    points between 0 and L.


### Version 0.50-18, 2006-12-18

-   A newly created package.

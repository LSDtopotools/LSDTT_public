The Land Surface Dynamics Chi package
=========================================

This is the code for performing chi analysis as described in this paper:

Mudd, S.M., M. Attal, D.T. Milodowski, S.W.D. Grieve and D.A. Valters (2014), A statistical framework to quantify spatial variation in channel gradients using the integral method of channel profile analysis, Journal of Geophysical Research-Earth Surface, 119 138-152, doi: 10.1002/2013JF002981 enhanced html.

[Click here to go to the open access version of the paper](http://onlinelibrary.wiley.com/enhanced/doi/10.1002/2013JF002981/)

If you have questions about this software please email simon.m.mudd _at_ ed.ac.uk

System requirements
======================================

This software works best in a linux environment. 

If you do not have a linux operating system, you can set one up in a [virtual machine](http://en.wikipedia.org/wiki/Virtual_machine).

Instructions for setting up a linux (Ubuntu) virtual machine on a windows operating sysem are [here.](http://www.geos.ed.ac.uk/~smudd/NMDM_Course/html/outside_edin.html)

You can also use [cygwin](https://www.cygwin.com/) to set up a unix-like environment on your windows computer. 
You will need to install the g++, gdb, and make utilities from the 'devel' menu. 

If you are using a fruit-based operating system, you can try [VirtualBox](https://www.virtualbox.org/). 

Compiling the code
=======================================

Go into the folder `chi_driver`.

Then use `make` to compile the code.

For example:

```
make -f chi_step1_write_junctions.make

make -f chi_step2_write_channel_file_driver.make

make -f chi_m_over_n_analysis.make

make -f chi_get_profiles.make
```

Preparing your data
===================================

At the moment the preferred data format is ENVI bil format because it retains georeferencing information. 
[Click here for details](http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/float.html)


Running the analysis
=====================================

You should read the following instructions. 
One difference between this version and the version in the documentation: the `make` and `driver` files are in the subfolder `chi_driver` rather than `driver_functions`.

[part1](http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/chi_profiles.html)

[part2](http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/chi_analysis_part2.html)

Note on constraining m/n
-----------------------------------

The Land Surface Dynamics Group at the University of Edinburgh has now run many hundreds of these analyses and have found that the collinearity test is more robust than single channel tests of the m/n ratio. 
In most cases, apart from landscapes with many hanging valleys, m/n should be determined with the collinearity test. 


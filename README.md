splitImageFract
=============

ETNA package attached to the paper:

_Andrea Azzarelli and Alessandro Buccini_       
**A variational model for glyph extraction, theoretical aspects and implementation**       
ETNA Volume ...

---

AUTHORS
-------

- Andrea Azzarelli,
       University of Cagliari, Italy
- Alessandro Buccini,
       University of Cagliari, Italy


> Email: andrea.azzarelli@unica.it, alessandro.buccini@unica.it

---

DESCRIPTION
-----------

This package provides one function capable of exctracting a glyph
embedded in an underlying background surface using ADMM and Fractional 
Laplacian. 

Release: 1.0, ...

Programming language: Matlab 25.1 (R2025a)

License: (see license.md)

---

DISCLAIMER
----------

ETNA editors have test run the code, but the reponsibility of the correctness
of the code is with the authors. ETNA is not responsible for any problems the 
code may cause; see the file License.md for details.

---

INSTALLATION
------------

Download and extract the splitImageFract package. 
To run the Test.m file the data files Circumference.mat, Triangle.mat, 
Writings.mat and DomusDeJanas.mat must be stored in the same folder or in the path. 
he files Test.m  contains the numerical examples present in the paper.

---

PACKAGE USE
------------

As long as all necessary files are in the working directory (see above), the
Test.m  demo may be executed without entry of additional information (i.e. the 
demos are already set-up).

---

PACKAGE STRUCTURE
-----------------

The following is a list with a brief description of the contents of the
splitImageFract package. For a more detailed description see the codePrimer.pdf
file in this package.


* License.md             : Markdown file containing license for package use. 
* README.md              : This file.
* splitImageFract.m	 :  the function implementing the splitting
* Test.m                 : A demo to showcase the use of the splitImageFract algorithm 
                           	on four examples.
* CreateData.m : A function to add rougnhess to the background
* Circumference.mat  :Synthetic   Dataset for the example 1
* Triangle.mat    : Synthetic Dataset for the example 2
* Writings.mat  : Synthetic Dataset for the example 3
* DomusDeJanas.mat  : Real Dataset for the example 4

In Synthetic Datasets are contained the variables
-G_true: exact solution
-S_true: Background smooth surface
-nu: parameter for CreateData.m
-sigma: parameter for CreateData.m
-mu: vector of parameters to be tested
-CLIM: for better visualization
-options: (optional) to change default parameters for splitImageFract

In the real Dataset are present the following parameters:
-D:  the surface to be splitted
-mu: vector of parameters to be tested

---


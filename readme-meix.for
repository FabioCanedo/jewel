		****************************************
		*           DISK TO ACCOMPANY          *
		*   COMPUTATION OF SPECIAL FUNCTIONS   *
		*                                      *
		*   Shanjie Zhang and Jianming Jin     *
		*                                      *
		*   Copyright 1996 by John Wiley &     *
		*              Sons, Inc.              *
		*                                      *
		****************************************

I. INTRODUCTION

     As stated in the preface of our book "Computation of Special 
Functions,"  the purpose of this book is to share with the reader  
a set of computer programs (130 in total) which we have developed 
during the past several years for computing a variety of  special  
mathematical functions.  For your convenience, we attach to the
book this diskette that contains all the computer programs  
listed or mentioned in the book. 

     In this diskette,  we place all the programs under directory 
SMF\PROGRAMS. In order to illustrate the use of these programs 
and facilitate your testing of the programs, we wrote a short 
simple main program for each program so that you can readily test 
them.    

     All the programs are written in FORTRAN-77 and tested on PCs
and workstations. Therefore, they should run on any computer with 
implementation of the FORTRAN-77 standard.  

     Although we have made a great effort to test these programs,  
we would not be surprised  to find some errors in them.  We would 
appreciate it if you can bring to our attention any errors you find.
You can do this by either writing us directly at the location
(e-mail: j-jin1@uiuc.edu) or writing to the publisher, whose address 
appears on the back cover of the book.  However, we must note that
all these programs are sold "as is," and we cannot guarantee to 
correct the errors reported by readers on any fixed schedule.

     All the programs and subroutines  contained in this diskette 
are copyrighted.   However,  we give permission to the reader who
purchases this book to incorporate any of these programs into his
or her programs provided that the copyright is acknowledged. 

     Regarding the specifics of the programs, we want to make the
following two points.

  1) All the programs are written in double precision.   Although
     the use of double precision is  necessary for some programs,
     especially for those based on series expansions,  it is not 
     necessary for all programs.  For example, the computation of
     of special functions based on polynomial approximations does
     not have to use double precision.  We chose to write all the
     programs using double precision in order to  avoid  possible 
     confusion  which  may  occur in using  these  programs.  If 
     necessary,  you can  convert  the programs into the single
     precision format easily.  However,  doing so for some 
     programs may lead to a lower accuracy.

  2) In the main programs that calculate a  sequence  of  special 
     functions, we usually set the maximum order or degree to 100
     or 250.  However, this is not a limit.  To compute functions  
     with a higher order or degree, all you need to do is simply 
     set the dimension of proper arrays higher.  


II. DISCLAIMER OF WARRANTY

     Although we have made a great effort to test and validate the 
computer programs, we make no warranties, express or implied, that 
these  programs  are  free  of  error,  or are consistent with any 
particular  standard  of  merchantability,  or that they will meet 
your requirements for any particular application.  They should not 
be relied on for  solving problems  whose incorrect solution could 
result in  injury to  a person or loss of property.  If you do use 
the programs in such a manner, it is at your own risk. The authors 
and publisher  disclaim all liability for  direct or consequential 
damages resulting from your use of the programs.


III. LIST OF PROGRAMS

(Please note that all file names of programs installed from the disk 
begin with an M, for example, MBERNOA.FOR)

BERNOA  Evaluate a sequence of Bernoulli numbers (method 1).

BERNOB  Evaluate a sequence of Bernoulli numbers (method 2).

EULERA  Evaluate a sequence of Euler numbers (method 1).

EULERB  Evaluate a sequence of Euler numbers (method 2). 

*****

OTHPL   Evaluate a sequence of orthogonal polynomials and their 
derivatives, including Chebyshev, Laguerre, and Hermite 
polynomials. 

LEGZO   Evaluate the nodes and weights for Gauss-Legendre quadrature.

LAGZO   Evaluate the nodes and weights for Gauss-Laguerre quadrature.

HERZO   Evaluate the nodes and weights for Gauss-Hermite quadrature.              

*****

GAMMA   Evaluate the gamma function.

LGAMA   Evaluate the gamma function or the logarithm of the gamma 
function.

CGAMA   Evaluate the gamma function with a complex argument.

BETA    Evaluate the beta function.

PSI     Evaluate the psi function.

CPSI    Evaluate the psi function with a complex argument.

INCOG   Evaluate the incomplete gamma function.

INCOB   Evaluate the incomplete beta function.

*****

LPN     Evaluate a sequence of Legendre polynomials and their 
derivatives with real arguments.

CLPN    Evaluate a sequence of Legendre polynomials and their 
derivatives with complex arguments.

LPNI    Evaluate a sequence of Legendre polynomials, their 
derivatives, and their integrals.

LQNA    Evaluate a sequence of Legendre functions of the second 
kind and their derivatives with restricted real arguments.

LQNB    Evaluate a sequence of Legendre functions of the second 
kind and their derivatives with nonrestricted real arguments.

CLQN    Evaluate a sequence of Legendre functions of the second 
kind and their derivatives with complex arguments.

LPMN    Evaluate a sequence of associated Legendre polynomials and 
their derivatives with real arguments.

CLPMN   Evaluate a sequence of associated Legendre polynomials and 
their derivatives with complex arguments.

LQMN    Evaluate a sequence of associated Legendre functions of the 
second kind and their derivatives with real arguments.

CLQMN   Evaluate a sequence of associated Legendre functions of the 
second kind and their derivatives with complex arguments.

LPMV    Evaluate associated Legendre functions of the first kind 
with an integer order and arbitrary non-negative degree.
 
*****

JY01A   Evaluate the zeroth- and first-order Bessel functions of the 
first and second kinds with real arguments using series and 
asymptotic expansions.

JY01B   Evaluate the zeroth- and first-order Bessel functions of the 
first and second kinds with real arguments using polynomial 
approximations.

JYNA    Evaluate a sequence of Bessel functions of the first and 
second kinds and their derivatives with integer orders and 
real arguments (method 1).

JYNB    Evaluate a sequence of Bessel functions of the first and 
second kinds and their derivatives with integer orders and 
real arguments (method 2).

CJY01   Evaluate the zeroth- and first-order Bessel functions of the 
first and second kinds and their derivatives with complex 
arguments.

CJYNA   Evaluate a sequence of Bessel functions of the first and 
second kinds and their derivatives with integer orders and 
complex arguments (method 1).

CJYNB   Evaluate a sequence of Bessel functions of the first and 
second kinds and their derivatives with integer orders and 
complex arguments (method 2).

JYV     Evaluate a sequence of Bessel functions of the first and 
second kinds and their derivatives with arbitrary real orders 
and real arguments.

CJYVA   Evaluate a sequence of Bessel functions of the first and 
second kinds and their derivatives with arbitrary real orders 
and complex arguments (method 1).

CJYVB   Evaluate a sequence of Bessel functions of the first and 
second kinds and their derivatives with arbitrary real orders 
and complex arguments (method 2).

CJK     Evaluate the coefficients for the asymptotic expansion of 
Bessel functions for large orders.

CJYLV   Evaluate Bessel functions of the first and second kinds and 
their derivatives with a large arbitrary real order and complex 
arguments.

JYZO    Evaluate the zeros of the Bessel functions of the first and 
second kinds and their derivatives.

JDZO    Evaluate the zeros of the Bessel functions of the first kind 
and their derivatives.

CYZO    Evaluate the complex zeros of the Bessel functions of the 
second kind of order zero and one.

LAMN    Evaluate a sequence of lambda functions with integer orders 
and their derivatives.

LAMV    Evaluate a sequence of lambda functions with arbitrary orders 
and their derivatives.

*****

IK01A   Evaluate the zeroth- and first-order modified Bessel 
functions of the first and second kinds with real arguments.

IK01B   Evaluate the zeroth- and first-order modified Bessel 
functions of the first and second kinds with real arguments.

IKNA    Evaluate a sequence of modified Bessel functions of the 
first and second kinds and their derivatives with integer 
orders and real arguments (method 1).

IKNB    Evaluate a sequence of modified Bessel functions of the 
first and second kinds and their derivatives with integer 
orders and real arguments (method 2).

CIK01   Evaluate the zeroth- and first-order modified Bessel 
functions of the first and second kinds and their derivatives 
with complex arguments.

CIKNA   Evaluate a sequence of modified Bessel functions of the 
first and second kinds and their derivatives with integer 
orders and complex arguments (method 1).

CIKNB   Evaluate a sequence of modified Bessel functions of the 
first and second kinds and their derivatives with integer 
orders and complex arguments (method 2).

IKV     Evaluate a sequence of modified Bessel functions of the 
first and second kinds and their derivatives with arbitrary 
real orders and real arguments.

CIKVA   Evaluate a sequence of modified Bessel functions of the 
first and second kinds and their derivatives with arbitrary 
real orders and complex arguments. 

CIKVB   Evaluate a sequence of modified Bessel functions of the 
first and second kinds and their derivatives with arbitrary 
real orders and complex arguments. 

CIKLV   Evaluate modified Bessel functions of the first and second 
kinds and their derivatives with a large arbitrary real order 
and complex arguments. 

CH12N   Evaluate a sequence of Hankel functions of the first and 
second kinds and their derivatives with integer orders and 
complex arguments.

*****

ITJYA   Evaluate the integral of Bessel functions J0(t) and Y0(t) 
from 0 to x using series and asymptotic expansions.

ITJYB   Evaluate the integral of Bessel functions J0(t) and Y0(t) 
from 0 to x using polynomial approximations.

ITTJYA  Evaluate the integral of [1-J0(t)]/t from 0 to x and Y0(t)/t   
from x to infinity using series and asymptotic expansions.

ITTJYB  Evaluate the integral of [1-J0(t)]/t from 0 to x and Y0(t)/t   
from x to infinity using polynomial approximations.

ITIKA   Evaluate the integral of modified Bessel functions I0(t) and   
K0(t) from 0 to x using series and asymptotic expansions.

ITIKB   Evaluate the integral of modified Bessel functions I0(t) and   
K0(t) from 0 to x using polynomial approximations.

ITTIKA  Evaluate the integral of [1-I0(t)]/t from 0 to x and K0(t) 
from x to infinity using series and asymptotic expansions.

ITTIKB  Evaluate the integral of [1-I0(t)]/t from 0 to x and K0(t) 
from x to infinity using polynomial approximations.

****

SPHJ    Evaluate a sequence of spherical Bessel functions of the 
first kind and their derivatives with integer orders and 
real arguments.

SPHY    Evaluate a sequence of spherical Bessel functions of the 
second kind and their derivatives with integer orders and 
real arguments.

CSPHJY  Evaluate a sequence of spherical Bessel functions of the 
first and second kinds and their derivatives with integer 
orders and complex arguments.

RCTJ    Evaluate a sequence of Riccati-Bessel functions and their 
derivatives of the first kind.

RCTY    Evaluate a sequence of Riccati-Bessel functions and their 
derivatives of the second kind.

SPHI    Evaluate a sequence of modified spherical Bessel functions 
of the first kind and their derivatives with integer orders 
and real arguments.

SPHK    Evaluate a sequence of modified spherical Bessel functions 
of the second kind and their derivatives with integer orders 
and real arguments.

CSPHIK  Evaluate a sequence of modified spherical Bessel functions 
of the first and second kinds and their derivatives with 
integer orders and complex arguments.

*****


KLVNA   Evaluate the Kelvin functions and their derivatives using 
series and asymptotic expansions.

KLVNB   Evaluate the Kelvin functions and their derivatives using 
polynomial approximations.

KLVNZO  Evaluate the zeros of the Kelvin functions and their 
derivatives.

***** 

AIRYA   Evaluate the Airy functions and their derivatives by means 
of Bessel functions.

AIRYB   Evaluate the Airy functions and their derivatives using the 
series and asymptotic expansions.

ITAIRY  Evaluate the integral of the Airy functions. 

AIRYZO  Evaluate the zeros of Airy functions and their derivatives.

***** 

STVH0   Evaluate the zeroth-order Struve function.

STVH1   Evaluate the first-order Struve function.

STVHV   Evaluate the Struve functions with an arbitrary order.

ITSH0   Evaluate the integral of Struve function H0(t) from 0 to x.

ITTH0   Evaluate the integral of H0(t)/t from x to infinity.

STVL0   Evaluate the zeroth-order modified Struve function.

STVL1   Evaluate the first-order modified Struve function.

STVLV   Evaluate the modified Struve function with an arbitrary 
order.

ITSL0   Evaluate the integral of modified Struve function L0(t) 
from 0 to x.

*****

HYGFX   Evaluate the hypergeometric function with real arguments.

HYGFZ   Evaluate the hypergeometric function with complex arguments.

*****

CHGM    Evaluate the confluent hypergeometric function M(a,b,x) with 
real arguments.

CCHG    Evaluate the confluent hypergeometric function M(a,b,z) with 
complex arguments.

CHGU    Evaluate the confluent hypergeometric function U(a,b,x) with 
real arguments.

*****

PBDV    Evaluate a sequence of parabolic cylinder functions Dv(x) and 
their derivatives.

PBVV    Evaluate a sequence of parabolic cylinder functions Vv(x) and 
their derivatives.

PBWA    Evaluate parabolic cylinder functions W(a,+/-x) and their 
derivatives.

CPBDN   Evaluate a sequence of parabolic cylinder functions Dn(z) and 
their derivatives for complex arguments.

***** 

CVA1    Evaluate a sequence of characteristic values for the Mathieu 
and modified Mathieu functions.

CVA2    Evaluate a specific characteristic value for the Mathieu 
and modified Mathieu functions.

FCOEF   Evaluate the expansion coefficients for the Mathieu and
modified Mathieu functions.

MTU0    Evaluate the Mathieu functions and their derivatives.

MTU12   Evaluate the modified Mathieu functions of the first and 
second kinds and their derivatives.

*****        

SEGV    Evaluate a sequence of characteristic values for spheroidal 
wave functions.

SDMN    Evaluate the expansion coefficients d_k^mn for spheroidal 
wave functions.

SCKA    Evaluate the expansion coefficients c_2k^mn for spheroidal 
wave functions (method 1).

SCKB    Evaluate the expansion coefficients c_2k^mn for spheroidal 
wave functions (method 2).

ASWFA   Evaluate the angular spheroidal wave functions of the first 
kind (method 1).

ASWFB   Evaluate the angular spheroidal wave functions of the first 
kind (method 2).

RSWFP   Evaluate the radial prolate spheroidal wave functions of the 
first and second kinds.

RSWFO   Evaluate the radial oblate spheroidal wave functions of the 
first and second kinds.

LPMNS   Evaluate a sequence of the associated Legendre functions of 
the first kind and their derivatives with real arguments  
for a given order.

LQMNS   Evaluate a sequence of the associated Legendre functions of 
the second kind and their derivatives with real arguments  
for a given order.

*****

ERROR   Evaluate the error function.

CERROR  Evaluate the error function with a complex argument.

*****

FCS     Evaluate the Fresnel Integrals.

FFK     Evaluate the modified Fresnel integrals.

CERZO   Evaluate the complex zeros of the error function.

FCSZO   Evaluate the complex zeros of the Fresnel Integrals.

***** 

CISIA   Evaluate the cosine and sine integrals using their series 
and asymptotic expansions.

CISIB   Evaluate the cosine and sine integrals using their rational 
approximations.

***** 

COMELP  Evaluate the complete elliptic integrals of the first and 
second kinds.

ELIT    Evaluate the incomplete elliptic integrals of the first and 
second kinds.

ELIT3   Evaluate the complete and incomplete elliptic integrals of 
the third kind.

JELP    Evaluate the Jacobian elliptic functions.

*****

E1XA    Evaluate the exponential integral E1(x) using its polynomial 
approximations.

E1XB    Evaluate the exponential integral E1(x) using its series and
continued fraction expressions.

E1Z     Evaluate the exponential integral E1(z) for complex arguments.

ENXA    Evaluate a sequence of exponential integrals En(x) (method 1).

ENXB    Evaluate a sequence of exponential integrals En(x) (method 2).

EIX     Evaluate the exponential integral Ei(x).





 
 

Code for magnetic field synthesis without singularities at the poles
====================================================================

This library provides routines to compute the magnetic flux vector based on
spherical harmonic coefficients. Adding the enclosing folder and its
subfolders to one’s path in Matlab will make all of the functions available
to the user.

All functions are implemented in Matlab. However, key functions are also
implemented in C++. The functions will execute significantly faster if the
CompileCLibraries.m function is run before using the library. A C++ compiler
must be installed and set up for the function CompileCLibraries.m to work.

The function DemonstrateMagneticCode.m demonstrates the basic features of how
the spherical harmonic coefficients are synthesized and is a good place to
start. The file Mathematical Function/spherHarmonicEval.m is the interface to
the main method for handling the spherical harmonic coefficients. The
different algorithms used to avoid singularities at the poles and a loss of
precision near the equator are documented therein. The function
spherHarmonicEval can also be used to synthesis gravitational coefficients
and coefficients for evaluating terrain heights.

Additionally, the functions Magnetism/geogHeading2Mag.m and Magnetism
magHeading2Geog.m can be used to convert between geographic and magnetic
headings to within the precision of the magnetic model.

If one does not wish to use the Matlab code, the C++ code that does not form
an interface to Matlab is generally broken out and kept in folders labeled
“Shared C++ Code”. Not all functions are implemented in C++. The key function 
spherHarmonicEval and the functions that it calls have C++ implementations.

The code has been tested under Matlab2013b on Windows 7 and Mac OS X 10.9.1.

LICENSE:

The source code is in the public domain and not licensed or under copyright.
The information and software may be used freely by the public. As required by
17 U.S.C. 403, third parties producing copyrighted works consisting
predominantly of the material produced by U.S. government agencies must provide
notice with such work(s) identifying the U.S. Government material incorporated
and stating that such material is not subject to copyright protection.

Derived works shall not identify themselves in a manner that implies an
endorsement by or an affiliation with the Naval Research Laboratory.

RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE SOFTWARE
AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL RESEARCH
LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT
IN THE USE OF THE SOFTWARE.

February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
E-Mail: david.crouse@nrl.navy.mil
In any published work or commercial product that uses the software in this
library, acknowledgement is appreciated.
(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
====================================================================
Files:
CompileCLibraries.m	A Matlab function that, when run, will compile all of
			the C++ implementations of the code. Due to the slow
			speed of the implementations done directly in Matlab,
			it is suggested that this function be run before
			using the library. The compiled code is placed in the
			folder Compiled_Code.
Constants.m		A class that collects the constants (other than
			magnetic field coefficients) that are used in the
			models. The class never has to be instantiated. the
			constants can all be access with the notation
			Constant.constantName.
DemonstrateMagneticCode.m A Matlab script that, when run, will evaluate the
			magnetic flux of the world magnetic model on a grid
			of points and then compute the offset of the
			direction of compass North from true North in degrees
			East of North (clockwise). The inclination angle of
			the magnetic flux vector is also plotted.
Compiled_Code 		(Folder) An empty folder in which the code compiled
			by the function CompileCLibraries is placed.
Container Classes	(Folder)
—ClusterSet.m		A class that allows sets of data to be easily
			accessed like a matrix.
—Shared C++ Code	(Folder)
—-ClusterSetCPP.hpp	A C++ class (implemented entirely in the header file)
			that overloads operators to provide a functionality
			similar to the ClusterSet class in Matlab.
CoordinateSystems	(Folder)
-calcSpherJacob.cpp 	A C++ implementation of calcSpherJacob.m.
-calcSpherJacob.m	Calculate the Jacobian of the spherical coordinate
			system at a point specified in Cartesian coordinates.
-Cart2Ellipse.m		Convert from Cartesian to ellipsoidal coordinates.
-Cart2Sphere.m		Convert from Cartesian to spherical coordinates.
-ellips2Cart.m		Convert from ellipsoidal to Cartesian coordinates.
-ellips2Sphere.m	Convert from ellipsoidal to spherical coordinates.
-getENUAxes.m		Get unit vectors defining the axes of the East-North
			-Up Coordinate system.
-getENUAxes.cpp		A C++ implementation of getENUAxes.m.
-rotateVector.m		Rotate a vector a specified angle about a given axis.
-spher2Cart.cpp		A C++ implementation of spher2Cart.m.
-spher2Cart.m		Convert from spherical to Cartesian coordinates.
-spher2Ellipse.m	Convert from spherical to ellipsoidal coordinates.
-Shared C++ Code	(Folder)
—-calcSpherJacobCPP.cpp A C++-only version of the function calcSpherJacob.
—-CoordFuncs.hpp	A C++ header for the functions in this folder.
—-getENUAxesCPP.cpp	A C++-only implementation of getENUAxes.
—-spher2CartCPP.cpp	A C++-only version of spher2Cart.
Gravity and Magnetism	(Folder)
-geogHeading2Mag.m	Convert from a direction given in radians East of
			geographic North to one in radians East of magnetic
			North.
-getWMM2010Coeffs.m	Load the coefficients for the 2010 World Magnetic
			Model.
-magHeading2Geog.m	Convert from a direction given in radians East of
			magnetic North to one in radians East of geographic
			North.
-data			(Folder)
—-WMM.COF		Coefficients for the 2010 World Magnetic Model.
Mathematical Functions	(Folder)
-NALegendreCosRat.cpp	A C++ implementation of NALegendreCosRat.m.
-NALegendreCosRat.m	A function to compute certain associated Legendre
			function ratios and derivatives for spherical
			harmonic synthesis.
-normHelmholtz.cpp	A C++ implementation of normHelmholtz.m
-normHelmholtz.m	A function to compute fully Normalized Helmholtz
			polynomials for spherical harmonic synthesis.
-spherHarmonicEval.m	A function to obtain potentials and gradients of
			potentials given spherical harmonic coefficients.
-spherHarmonicEvalCPPInt.cpp A function used internally to
			spherHarmonicEval.m that interfaces it with a faster
			C++ implementation of the function.
-Shared C++ Code	(Folder)
—-mathFuncs.hpp		A C++ header files for the functions implemented in
			this folder.
—-NALegendreCosRatCPP.cpp A C++-only version of the features of
			NALegendreCosRat.m.
—-normHelmholtzCPP.cpp	A C++-only version of normHelmholtz.m.
—-spherHarmonicEvalCPP.cpp A C++-only version of spherHarmonicEval.m.
Misc			(Folder)
-MexValidation.h	A set of functions for simplifying the interface
			between Matlab and C/ C++ code.
-systemNumberOfBits.m	A function that says whether the system is 32 or 64
			bit based on how Matlab named its compiled code.
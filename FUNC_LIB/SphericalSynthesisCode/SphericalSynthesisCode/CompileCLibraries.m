%%COMPILECLIBRARIES Compile all of the functions implemented in C++ that
%                   are to be called from Matlab. If a C or C++ mex
%                   function has the same name as a file implemented in
%                   Matlab and is compiled, then Matlab will execute the
%                   compiled mex function rather than the (Usually slower)
%                   Matlab code.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

function CompileCLibraries()

%Get the path to this file. it should be root level in the library.
ScriptPath=mfilename('fullpath');
ScriptFolder = fileparts(ScriptPath);

%Save the old working directory and switch to the location of this file.
curDir=pwd;
cd(ScriptFolder)

%Compile general coordinate system code.
mex('-v','-largeArrayDims','-outdir','./Compiled_Code/','-I./Misc/','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/spher2Cart.cpp','./Coordinate Systems/Shared C++ Code/spher2CartCPP.cpp');
mex('-v','-largeArrayDims','-outdir','./Compiled_Code/','-I./Misc/','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/calcSpherJacob.cpp','./Coordinate Systems/Shared C++ Code/calcSpherJacobCPP.cpp','./Coordinate Systems/Shared C++ Code/spher2CartCPP.cpp');
mex('-v','-largeArrayDims','-outdir','./Compiled_Code/','-I./Misc/','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/getENUAxes.cpp','./Coordinate Systems/Shared C++ Code/getENUAxesCPP.cpp');

%Compile the magnetic code.
mex('-v','-largeArrayDims','-outdir','./Compiled_Code/','-I./Misc/','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Mathematical Functions/NALegendreCosRat.cpp','./Mathematical Functions/Shared C++ Code/NALegendreCosRatCPP.cpp');
mex('-v','-largeArrayDims','-outdir','./Compiled_Code/','-I./Misc/','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Mathematical Functions/normHelmholtz.cpp','./Mathematical Functions/Shared C++ Code/normHelmholtzCPP.cpp');
mex('-v','-largeArrayDims','-outdir','./Compiled_Code/','-I./Misc/','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','-I./Coordinate Systems/Shared C++ Code/','./Mathematical Functions/spherHarmonicEvalCPPInt.cpp','./Mathematical Functions/Shared C++ Code/spherHarmonicEvalCPP.cpp','./Mathematical Functions/Shared C++ Code/NALegendreCosRatCPP.cpp','./Mathematical Functions/Shared C++ Code/normHelmholtzCPP.cpp','./Coordinate Systems/Shared C++ Code/spher2CartCPP.cpp','./Coordinate Systems/Shared C++ Code/calcSpherJacobCPP.cpp');

%Restore the old working directory.
cd(curDir);
end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
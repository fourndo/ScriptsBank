classdef Constants
%%CONSTANTS Static methods for obtaining standard physical constants.
%
%Constants can be accessed using Constants.constantName without
%instantiating this class.
%
%The WGS84 properties are from 
%Department of Defense, "Department of Defense world geodetic system 1984:
%Its definition and relationships with local geodetic systems," National
%Imagery and Mapping Agency, Tech. Rep. NIMA TR8350.2, Jun. 2004, third
%Edition, Amendment 2. [Online]. Available:
%http://earth- info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
%
%The ellipsoid properties for the EGM2008 model are taken from 
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/README_FIRST.pdf
%where the value of GM and the semi-major axis of the reference ellipsoid
%differ from the WGS84 values. Note that these values are the same for the
%EGM96 gravitational model.
%
%The value of the reference sphere radius in the World Magnetic Model for
%the year 2010 is given in
%S. Maus, S. McLean, M. Nair, and C. Rollins, "The US/UK
%world magnetic model for 2010-2015," National Oceanographic and
%Atmospheric Organization, Tech. Rep. NESDIS/NGDC, 2010. [Online].
%Available: http://www.ngdc.noaa.gov/geomag/WMM/
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

properties (Constant)
    %WGS84 Properties
    %GM is the universal gravitation constant times the mass of the Earth.
    WGS84GMWithAtmosphere=3.986004418*10^(14);%m^3/s^2
    WGS84GMWithoutAtmosphere=3.9860009*10^(14);%m^3/s^2
    WGS84EarthRotationRate=7292115.0*10^(-11);%radians per second.
    %The following 3 parameters should be consistent with values in the
    %IAU's SOFA library.
    WGS84SemiMajorAxis=6378137.00;%m
    WGS84InverseFlattening=298.257223563;%Unitless
    WGS84Flattening=1/298.257223563;%Unitless
    
    %The EGM2008 model; the same values are used for the EGM96 model. These
    %values are needed for using the spherical harmonic coefficients in
    %the models. They are also the defining parameters of the reference
    %ellipsoid used in the model
    EGM2008GM=3986004.415*10^8;%m^3/s^2
    EGM2008SemiMajorAxis=6378136.3%m
    EGM2008EarthRotationRate=7292115*10^(-11);%rad/s
    EGM2008C20Bar=-484.1654767*10^(-6);%Unitless, tide-free, defines the
    %reference ellipsoid
    
    %The radius of the reference sphere used in the WMM2010.
    WMM2010SphereRad=6371200;%meters
end
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

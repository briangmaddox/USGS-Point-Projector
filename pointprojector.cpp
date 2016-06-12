
//pointprojector projects some points from one specified projection to another


#include "ProjectionLib/UTMProjection.h"
#include "ProjectionLib/StatePlaneProjection.h"
#include "ProjectionLib/AlbersConicProjection.h"
#include "ProjectionLib/LambertConformalConicProjection.h"
#include "ProjectionLib/MercatorProjection.h"
#include "ProjectionLib/PolarStereographicProjection.h"
#include "ProjectionLib/PolyconicProjection.h"
#include "ProjectionLib/EquidistantConicProjection.h"
#include "ProjectionLib/TransverseMercatorProjection.h"
#include "ProjectionLib/StereographicProjection.h"
#include "ProjectionLib/LambertAzimuthalProjection.h"
#include "ProjectionLib/AzimuthalEquidistantProjection.h"
#include "ProjectionLib/GnomonicProjection.h"
#include "ProjectionLib/OrthographicProjection.h"
#include "ProjectionLib/SinusoidalProjection.h"
#include "ProjectionLib/EquirectangularProjection.h"
#include "ProjectionLib/MillerCylindricalProjection.h"
#include "ProjectionLib/VanDerGrintenProjection.h"
#include "ProjectionLib/RobinsonProjection.h"
#include "ProjectionLib/AlaskaConformalProjection.h"
#include "ProjectionLib/HotineObliqueMercatorProjection.h"
#include "ProjectionMesh/ProjectionMesh.h"
#include <getopt.h>
#include "math.h"
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>

using namespace ProjLib;

Projection * SetProjection(std::string parameterfile) throw();
DATUM GetDatum(std::string indatum) throw();
UNIT  GetUnit(std::string inunit) throw();
void getMinMax(std::vector<double>& array, double& min, double& max)
  throw();
double ConvertToDMS(double degrees) throw();

int main(int argc, char ** argv)
{
  Projection * to, * from;              //projections
  PmeshLib::ProjectionMesh pmesh; 
  std::string infile, inparam, outparam, 
    outfile("pout.txt");                //input/output files
  bool usepmesh = false;                //use the pmesh
  int pmeshsize = 4;                    //default pmesh size
  int myopt;                            //for options
  int interpolator = 7;                 //which interpolator for pmesh
  long int counter;
  double x, y;                          //for points
  long int numpoints;                   //numberofpoints
  std::ifstream in;
  std::ofstream out;
  std::ofstream compare;
  double min_x, min_y, max_x, max_y;
  double tempx, tempy;
  std::vector<double> xarr, yarr;
  double errorx, errory;
  std::string comparefile("compare.txt");

  try
    {
      //first parse the command args
       //parse the arguments
     while (( myopt = getopt(argc, argv, "c:n:l:p:f:i:o:?")) != -1)
       switch(myopt)
         {
         case 'c':
           {
             if (optarg)
               comparefile = optarg;
             break;
           }
         case 'n':
           {
             if (optarg)
               pmeshsize = atoi(optarg);
             break;
           }
         case 'p':
           {
             if (optarg)
               {
                 interpolator = atoi(optarg);
                 usepmesh = true;
               }
             break;
           }
         case 'f': //input file switch
           {
             if (optarg) // only change if they entered something
               infile = std::string(optarg);
             break;
           }
         case 'i': //parameter file for input
           {
             if (optarg)
               inparam = std::string(optarg);
             break;
           }
         case 'o':
           {
             if (optarg)
               outparam = std::string(optarg);
             break;
           }
         case '?':
         default: // help
           {
             std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
             std::cout << "where options are: " << std::endl;
             std::cout << "   -f input file to reproject" << std::endl;
             std::cout << "   -i input parameter file" << std::endl;
             std::cout << "   -? this help screen"   
                       << "(Which you obviously already knew about)"
                       << std::endl;
             return 0;;
           }
         }

     from = SetProjection(inparam);
     to = SetProjection(outparam);
     in.open(infile.c_str());
     
     in >> numpoints;
     
     xarr.resize(numpoints);
     yarr.resize(numpoints);
     
     for (counter = 0; counter < numpoints; counter++)
       in >> xarr[counter] >> yarr[counter];
     
     in.close();
     
     if (usepmesh)
       {
         // getMinMax(xarr, min_x, max_x);
         // getMinMax(yarr, min_y, max_y);
         min_x = 314850.0;
         max_x = min_x + 6224;
         max_y =  4291485.0;
         min_y = max_y - 7673.0;

          //setup the foward mesh
         pmesh.setSourceMeshBounds(min_x, min_y, max_x, max_y);
         pmesh.setMeshSize(pmeshsize,pmeshsize); 
    
         pmesh.setInterpolator(interpolator);
         //now project all the points onto the mesh
         pmesh.calculateMesh((*from), (*to));
     
       }
     
     out.open(outfile.c_str());
     compare.open(comparefile.c_str());

     for (counter = 0; counter < numpoints; counter++)
       {
         if (usepmesh)
           {
             tempx = xarr[counter];
             tempy = yarr[counter];
             from->projectToGeo(xarr[counter], yarr[counter], tempy, 
                                tempx);
             to->projectFromGeo(tempy, tempx, tempx, tempy);

             x = xarr[counter];
             y = yarr[counter];
             pmesh.projectPoint(x, y);
             
             errorx = tempx - x;
             errory = tempy - y;

             if (to->getUnit() == US_FEET)
               {
                 
                 errorx = errorx/3.28083988; //convert to meters
                 errory = errory/3.28083988; //convert to meters
               }
                 if (errorx < 0)
                   errorx = -errorx;
                 if (errory < 0)
                   errory = -errory;
                 
                 compare << tempx << "\t" << tempy
                         << "\t" << x << "\t" << y
                         << "\t" << errorx 
                         << "\t" << errory << "\t";
                 errorx = sqrt(errorx*errorx + errory*errory);
                 compare << errorx << std::endl; //distance
           }
         else
           {
             from->projectToGeo(xarr[counter], yarr[counter], y, x);
             to->projectFromGeo(y, x, x, y);
           }
         
         out << x << " " << y << std::endl;
       }
     out.close();
     compare.close();
     delete from;
     delete to;
     return 0;
    }
  catch(...)
    {
      delete from;
      delete to;
      return 0;
    }
}



Projection * SetProjection(std::string parameterfile) throw()
{
   
  std::ifstream in;
  std::string projtype, sdatum, sunit;
  double StdParallel1 = 0.0;
  double StdParallel2 = 0.0;
  double NatOriginLong = 0.0;
  double NatOriginLat = 0.0;
  double FalseEasting = 0.0;
  double FalseNorthing = 0.0;
  double FalseOriginLong = 0.0;
  double FalseOriginLat = 0.0;
  double FalseOriginEasting = 0.0;
  double FalseOriginNorthing = 0.0;
  double CenterLong = 0.0;
  double CenterLat = 0.0;
  double CenterEasting = 0.0;
  double CenterNorthing = 0.0;
  double ScaleAtNatOrigin = 1.0;
  double AzimuthAngle = 0.0;
  double StraightVertPoleLong = 0.0;
  int zone = 0;
  Projection * proj;
  
  try
    {
      in.open(parameterfile.c_str());
      
      
      in >> projtype;
      
      if (projtype == std::string("UTM"))
        {
          in >> zone;
          in >> sdatum;
          in >> sunit;
          
          
          if (!(proj = new (std::nothrow) UTMProjection(zone, 
                                                        GetDatum(sdatum), 
                                                        GetUnit(sunit))))
            throw std::bad_alloc();
        }
      
      if (projtype == std::string("SPCS"))
        {
          in >> zone;
          in >> sdatum;
          in >> sunit;
          StatePlaneProjection::setNAD83ParameterFilename
            (std::string("./nad83sp"));
          //create the state plane projection
          //here I am assuming that the zone will be 
          //a correct state plane zone
          if (!(proj = new (std::nothrow) StatePlaneProjection
                (zone, GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      
      if (projtype == std::string("ALBERS"))
        {
          in >> sdatum;
          in >> sunit;
          in >> StdParallel1
             >> StdParallel2
             >> CenterLong
             >> NatOriginLat
             >> FalseEasting
             >> FalseNorthing;
          
          if(!(proj =  new (std::nothrow) AlbersConicProjection
               ( ConvertToDMS(StdParallel1), 
                 ConvertToDMS(StdParallel2),
                 0.0, 0.0, ConvertToDMS(CenterLong),
                 ConvertToDMS(NatOriginLat), 
                 FalseEasting, FalseNorthing, GetDatum(sdatum), 
                 GetUnit(sunit))))
            throw std::bad_alloc();
        
          
        }
      if (projtype == std::string("AZMEQD"))
        {
          in >> sdatum >> sunit >> CenterLong >> CenterLat
             >> FalseEasting >> FalseNorthing;
          

          if(!(proj = new (std::nothrow) AzimuthalEquidistantProjection
               ( ConvertToDMS(CenterLong),
                 ConvertToDMS(CenterLat),
                 FalseEasting, FalseNorthing, 0.0, GetDatum(sdatum), 
                 GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("GNOMON"))
        {
          in >> sdatum >> sunit >> CenterLong >> CenterLat
             >> FalseEasting >> FalseNorthing;

          if(!(proj = new(std::nothrow) GnomonicProjection
               ( ConvertToDMS(CenterLong),
                 ConvertToDMS(CenterLat),
                 FalseEasting,FalseNorthing, 
                 0.0, GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("LAMAZ"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> CenterLat >> FalseEasting >> FalseNorthing;
          
          if(!(proj = new(std::nothrow) LambertAzimuthalProjection
               ( ConvertToDMS(CenterLong), ConvertToDMS(CenterLat),
                 FalseEasting, FalseNorthing, 0.0, 
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("ORTHO"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> CenterLat >> FalseEasting >> FalseNorthing;
          
          if(!(proj = new(std::nothrow) OrthographicProjection
               ( ConvertToDMS(CenterLong),
                 ConvertToDMS(CenterLat),
                 FalseEasting,FalseNorthing, 0.0, 
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("STEREO"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> CenterLat >> FalseEasting >> FalseNorthing;
          
          if(!(proj = new(std::nothrow) StereographicProjection
               (ConvertToDMS(CenterLong),
                ConvertToDMS(CenterLat),
                FalseEasting, FalseNorthing,
                0.0, GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("MILLER"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> FalseEasting >> FalseNorthing;
          
          if(!(proj = new (std::nothrow) MillerCylindricalProjection
               ( 0.0, ConvertToDMS(CenterLong),
                 FalseEasting, FalseNorthing,
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("ROBIN"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> FalseEasting >> FalseNorthing;

          if(!(proj = new(std::nothrow) RobinsonProjection
               ( 0.0, ConvertToDMS(CenterLong),
                 FalseEasting, FalseNorthing, 
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("SNSOID"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> FalseEasting >> FalseNorthing;
          if(!(proj = new(std::nothrow) SinusoidalProjection
               ( 0.0, ConvertToDMS(CenterLong),
                 FalseEasting,
                 FalseNorthing, GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("EQUIDC"))
        {
          in >> sdatum >> sunit >> StdParallel1
             >> StdParallel2;
          
          if ( StdParallel1 == StdParallel2 )
            {
              in >> CenterLat >> CenterLong >> NatOriginLat
                 >> FalseEasting >> FalseNorthing;
              if(!(proj =  new(std::nothrow) EquidistantConicProjection
                   ( ConvertToDMS(CenterLat), 0.0, 0.0,
                     ConvertToDMS(CenterLong),
                     ConvertToDMS(NatOriginLat),
                     FalseEasting,
                     FalseNorthing,
                     GetDatum(sdatum), GetUnit(sunit))))
                throw std::bad_alloc();
              
            }
          else
            {      
              in >> CenterLong >> NatOriginLat >> FalseEasting
                 >> FalseNorthing;
              if(!(proj = new(std::nothrow) EquidistantConicProjection
                   ( ConvertToDMS(StdParallel1),
                     ConvertToDMS(StdParallel2),
                     0.0, 0.0,
                     ConvertToDMS(CenterLong),
                     ConvertToDMS(NatOriginLat),
                     FalseEasting,
                     FalseNorthing,
                     GetDatum(sdatum), GetUnit(sunit))))
                throw std::bad_alloc();
            }
        }
      if (projtype == std::string("EQRECT"))
        {
          in >> sdatum >> sunit >> CenterLat >> CenterLong
             >> FalseEasting >> FalseNorthing;
          
          if(!(proj =  new(std::nothrow) EquirectangularProjection
               ( ConvertToDMS(CenterLat), 0.0, 
                 ConvertToDMS(CenterLong),
                 FalseEasting,FalseNorthing, 
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("HOM"))
        {
          in >> sdatum >> sunit >> ScaleAtNatOrigin
             >> AzimuthAngle >> CenterLong
             >> CenterLat >> FalseEasting 
             >> FalseNorthing;
          
          if(!(proj = new (std::nothrow) HotineObliqueMercatorProjection
               ( ScaleAtNatOrigin, AzimuthAngle,
                 0.0, 0.0, ConvertToDMS(CenterLong),
                 ConvertToDMS(CenterLat), FalseEasting,
                 FalseNorthing, GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("LAMCC"))
        {
          in >> sdatum >> sunit >> StdParallel1 >> StdParallel2
             >> NatOriginLong >> FalseOriginLat
             >> FalseEasting >> FalseNorthing;
          
          if(!(proj = new(std::nothrow) LambertConformalConicProjection
               ( ConvertToDMS(StdParallel1), ConvertToDMS(StdParallel2),
                 0.0, 0.0, ConvertToDMS(NatOriginLong),
                 ConvertToDMS(FalseOriginLat),
                 FalseEasting,FalseNorthing, GetDatum(sdatum), 
                 GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("MERCAT"))
        {
          in >> sdatum >> sunit >> NatOriginLong
             >> NatOriginLat >> CenterEasting 
             >> CenterNorthing;
          
          if(!(proj = new(std::nothrow) MercatorProjection
               ( 0.0, 0.0, ConvertToDMS(NatOriginLong),
                 ConvertToDMS(NatOriginLat),
                 CenterEasting, CenterNorthing, 
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("POLYC"))
        {
          
          in >> sdatum >> sunit >> CenterLong
             >> CenterLat >> FalseEasting
             >> FalseNorthing;
          
          if(!(proj = new(std::nothrow) PolyconicProjection
               ( 0.0, 0.0, ConvertToDMS(CenterLong),
                 ConvertToDMS(CenterLat), FalseEasting,
                 FalseNorthing, GetDatum(sdatum), 
                 GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("PS"))
        {
          in >> sdatum >> sunit >> StraightVertPoleLong
             >> NatOriginLat >> FalseEasting >> FalseNorthing;
          
          if(!(proj = new(std::nothrow) PolarStereographicProjection
               (ConvertToDMS(StraightVertPoleLong),
                ConvertToDMS(NatOriginLat), 0.0, 0.0,
                FalseEasting,FalseNorthing, 
                GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("ALASKA"))
        {
          in >> sdatum >> sunit >> FalseEasting 
             >> FalseNorthing;
          
          if(!(proj =  new(std::nothrow) AlaskaConformalProjection
                   (0.0, 0.0, FalseEasting, FalseNorthing,
                    GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("TM"))
        {
          in >> sdatum >> sunit >> ScaleAtNatOrigin 
             >> CenterLong >> NatOriginLat
             >> FalseEasting >> FalseNorthing;
          
          if(!(proj = new(std::nothrow) TransverseMercatorProjection
               (ScaleAtNatOrigin, 0.0, 0.0,
                ConvertToDMS(CenterLong),
                ConvertToDMS(NatOriginLat), FalseEasting,
                FalseNorthing, 
                GetDatum(sdatum), GetUnit(sunit))))
              throw std::bad_alloc();
        }
      if (projtype == std::string("VGRINT"))
        {
          in >> sdatum >> sunit >> CenterLat
             >> CenterLong >> FalseEasting 
             >> FalseNorthing;
          if(!(proj = new(std::nothrow) VanDerGrintenProjection
               ( ConvertToDMS(CenterLat), 0.0, 
                 ConvertToDMS(CenterLong),
                 FalseEasting,FalseNorthing,
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      in.close();
      return proj;
    }
  catch(...)
    {
      in.close();
      return NULL;
    }
}
      
      
DATUM GetDatum(std::string indatum) throw()
{
   
   if (indatum == std::string("ADINDAN"))
     return ADINDAN;
   
   if (indatum == std::string("ARC1950"))
     return ARC1950;
   
   if (indatum == std::string("ARC1960"))
     return ARC1960;
   
   if (indatum == std::string("AUSTRALIAN_GEODETIC_1966"))
     return AUSTRALIAN_GEODETIC_1966;
   
   if (indatum == std::string("AUSTRALIAN_GEODETIC_1984"))
     return AUSTRALIAN_GEODETIC_1984;
     
   if (indatum == std::string("CAPE"))
     return CAPE;
   
   if (indatum == std::string("EUROPEAN_DATUM_1950"))
     return EUROPEAN_DATUM_1950;
   
   if (indatum == std::string("HU_TZU_SHAN"))
     return HU_TZU_SHAN;
   
   if (indatum == std::string("INDIAN"))
     return INDIAN;

   if (indatum == std::string("NAD27"))
     return NAD27;
     
   if (indatum == std::string("NAD83"))
     return NAD83;

   if (indatum == std::string("ORDNANCE_SURVEY_1936"))
     return ORDNANCE_SURVEY_1936;

   if (indatum == std::string("PULKOVO_1942"))
     return PULKOVO_1942;

   if (indatum == std::string("PROVISIONAL_S_AMERICAN_1956"))
     return PROVISIONAL_S_AMERICAN_1956;

   if (indatum == std::string("TOKYO"))
     return TOKYO;
   
   if (indatum == std::string("WGS_72"))
     return WGS_72;
   
   if (indatum == std::string("WGS_84"))
     return WGS_84;
   
   return UNKNOWN_DATUM;
}
  
UNIT  GetUnit(std::string inunit) throw()
{
  if (inunit == std::string("METERS"))
    return METERS;
  
  if ((inunit == std::string("FEET")) || (inunit == std::string("US_FEET")))
    return US_FEET;

  return UNKNOWN_UNIT;

}

// **************************************************************************
double ConvertToDMS(double degrees) throw()
{
  double temp;
  int deg, min;
  
  int sign = 1;
  
  temp = degrees;
  if (degrees < 0.0)
  {
    sign = -1;
    temp = -temp;
  }
  //get the degrees
  deg = static_cast<int>(temp);
  temp -= deg;
  temp *= 60.0; //get minutes
  min = static_cast<int>(temp);
  temp -= min; //get seconds
  temp *= 60;

  temp = deg * 1000000 + min * 1000 + temp;
  return temp*sign;
}






// ************************************************************************
void getMinMax(std::vector<double>& array, double& min, double& max)
  throw()
{
  std::sort(array.begin(), array.end());

  max = array[array.size() - 1];
  min = array[0];
}

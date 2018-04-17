/** Derived from DTGeometryAnalyzer by Nicola Amapane
 *
 *  \author M. Maggi - INFN Bari
 *  \author A. Purohit - SINP Kolkata updated the code and changed the output format to HTML to make the output more user friendly.
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include <memory>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <set>

class ME0GeometryAnalyzer : public edm::one::EDAnalyzer<>
{
public: 
  ME0GeometryAnalyzer( const edm::ParameterSet& pset);

  ~ME0GeometryAnalyzer();

  void beginJob() override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const&) override;
  void endJob() override {}
  
private:
  const std::string& myName() { return myName_;}

  const int dashedLineWidth_;
  const std::string dashedLine_;
  const std::string myName_;
  std::ofstream ofos;
};

using namespace std;

ME0GeometryAnalyzer::ME0GeometryAnalyzer( const edm::ParameterSet& /*iConfig*/ )
  : dashedLineWidth_(104), dashedLine_( string(dashedLineWidth_, '-') ), 
    myName_( "ME0GeometryAnalyzer" ) 
{ 
  ofos.open("MytestOutput.html"); 
  ofos <<"<!DOCTYPE html>"<< endl;
  ofos <<"<html>"<< endl;
  ofos << "<head>" << endl << "<style>" << endl << "   table, th, td {" << endl << " border: 1px solid black;" << endl << " border-collapse: collapse;"<< endl << " }" << endl <\
    < " th {" << endl << "  padding: 5x;" << endl << "   text-align: center; ""   text-align: center; " << endl << "color:red;" << endl << "font-size:160%;" << endl <<  " td { " <<\
    endl << "  padding: 5px;" << endl <<"   text-align: center; "<< endl << " } " << endl << "</style>" << endl;
  ofos << "<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js\"></script>" << endl;
  ofos << "<script>" << endl
       << " $(document).ready(function(){" << endl
       << " $(\"#hide\").click(function(){" << endl
       << "  $(\"table\").hide();" << endl
       << "});" << endl
       << "     $(\"#show\").click(function(){" << endl
       << "  $(\"table\").show();" << endl
       << " });" << endl
       << " });" << endl
       << "</script>" << endl;
  ofos << "</head>" <<endl;
  ofos <<"<body>"<< endl;
  ofos <<"<p1> ======================== Opening output file </p1>"<< endl;

  // ofos <<"======================== Opening output file"<< endl;
}

ME0GeometryAnalyzer::~ME0GeometryAnalyzer() 
{
  ofos.close();
  ofos <<"======================== Closing output file"<< endl;
}

void
ME0GeometryAnalyzer::analyze( const edm::Event& /*iEvent*/, const edm::EventSetup& iSetup )
{
  edm::ESHandle<ME0Geometry> pDD;
  iSetup.get<MuonGeometryRecord>().get(pDD);     
  
  ofos << myName() << ": Analyzer..." << endl;
  ofos << "start " << dashedLine_ << endl;

  ofos << " Geometry node for ME0Geom is  " << &(*pDD) << endl;   
  ofos << " detTypes       \t"              <<pDD->detTypes().size() << endl;
  ofos << " GeomDetUnit       \t"           <<pDD->detUnits().size() << endl;
  ofos << " GeomDet           \t"           <<pDD->dets().size() << endl;
  ofos << " GeomDetUnit DetIds\t"           <<pDD->detUnitIds().size() << endl;
  ofos << " eta partitions \t"              <<pDD->etaPartitions().size() << endl;

  // checking uniqueness of roll detIds 
  bool flagNonUniqueRollID = false;
  bool flagNonUniqueRollRawID = false;
  for (auto roll1 : pDD->etaPartitions()){
    for (auto roll2 : pDD->etaPartitions()){
      if (roll1 != roll2){
	if (roll1->id() == roll2->id()) flagNonUniqueRollID = true;
	if (roll1->id().rawId() == roll2->id().rawId()) flagNonUniqueRollRawID = true;
      }
    }
  }
  // checking the number of strips and pads
  if (flagNonUniqueRollID or flagNonUniqueRollRawID)
    ofos << " -- WARNING: non unique roll Ids!!!" << endl;

  ofos << myName() << ": Begin iteration over geometry..." << endl;
  ofos << "iter " << dashedLine_ << endl;
  
  ofos << myName() << "Begin ME0Geometry TEST" << endl;
  ofos<<"<table style=\"width:100%\">" << endl << "<tr>" << endl;
  ofos << "<th> GEM Super Chamber Id </th>" << endl << "<th> Chamber </th>" << endl << "<th> Roll </th>" << endl << "<th> r (bottom) in cm </th>" << endl << "<th> W in cm </th>\
" << endl << "<th> h in cm </th>" << endl << "<th> cStrip1 phi </th>"<< endl << "<th> cStripN phi </th>"<< endl << "<th> dphi </th>"<< endl << "</tr>" << endl;

  /*
   * possible checklist for an eta partition:
   *   base_bottom, base_top, height, strips, pads
   *   cx, cy, cz, ceta, cphi
   *   tx, ty, tz, teta, tphi
   *   bx, by, bz, beta, bphi
   *   pitch center, pitch bottom, pitch top
   *   deta, dphi
   *   gap thicess
   *   sum of all dx + gap = chamber height
   */      
	    
  for (auto roll : pDD->etaPartitions()){
    ME0DetId rId(roll->id());
    ofos<<"\tME0EtaPartition , ME0DetId = " << rId.rawId() << ", " << rId << endl;
	      
    const BoundPlane& bSurface(roll->surface());
    const StripTopology* topology(&(roll->specificTopology()));
    
    // base_bottom, base_top, height, strips, pads (all half length)
    auto& parameters(roll->specs()->parameters());
    float bottomEdge(parameters[0]);
    float topEdge(parameters[1]);
    float height(parameters[2]);
    //    float nStrips(parameters[3]);
    //float nPads(parameters[4]);
    
    LocalPoint  lCentre( 0., 0., 0. );
    GlobalPoint gCentre(bSurface.toGlobal(lCentre));
    
    LocalPoint  lTop( 0., height, 0.);
    GlobalPoint gTop(bSurface.toGlobal(lTop));
    
    LocalPoint  lBottom( 0., -height, 0.);
    GlobalPoint gBottom(bSurface.toGlobal(lBottom));
    
    //   gx, gy, gz, geta, gphi (center)
    double cx(gCentre.x());
    double cy(gCentre.y());
    double cz(gCentre.z());
    double ceta(gCentre.eta());
    int cphi(static_cast<int>(gCentre.phi().degrees()));
    if (cphi < 0) cphi += 360;
    
    double tx(gTop.x());
    double ty(gTop.y());
    double tz(gTop.z());
    double teta(gTop.eta());
    int tphi(static_cast<int>(gTop.phi().degrees()));
    if (tphi < 0) tphi += 360;
	      
    double bx(gBottom.x());
    double by(gBottom.y());
    double bz(gBottom.z());
    double beta(gBottom.eta());
    int bphi(static_cast<int>(gBottom.phi().degrees()));
    if (bphi < 0) bphi += 360;
    
    /*
    // pitch bottom, pitch top, pitch centre
    float pitch(roll->pitch());
    float topPitch(roll->localPitch(lTop));
    float bottomPitch(roll->localPitch(lBottom));
    */
    // Type - should be ME0 Somethng
    string type(roll->type().name());
    
    // print info about edges
    LocalPoint lEdge1(topology->localPosition(0.));
    /*
    LocalPoint lEdgeN(topology->localPosition((float)nStrips));
    
    double cstrip1(roll->toGlobal(lEdge1).phi().degrees());
    double cstripN(roll->toGlobal(lEdgeN).phi().degrees());
    double dphi(cstripN - cstrip1);
    if (dphi < 0.) dphi += 360.;
    */
    double deta(abs(beta - teta));
    const bool printDetails(true);
    if (printDetails)
      ofos << "\t\tType: " << type << endl
	   << "\t\tDimensions[cm]: b = " << bottomEdge*2 << ", B = " << topEdge*2 << ", H  = " << height*2 << endl
	//	   << "    \tnStrips = " << nStrips << ", nPads =  " << nPads << endl
	   << "\t\ttop(x,y,z)[cm] = (" << tx << ", " << ty << ", " << tz << "), top (eta,phi) = (" << teta << ", " << tphi << ")" << endl
	   << "\t\tcenter(x,y,z) = (" << cx << ", " << cy << ", " << cz << "), center(eta,phi) = (" << ceta << ", " << cphi << ")" << endl
	   << "\t\tbottom(x,y,z) = (" << bx << ", " << by << ", " << bz << "), bottom(eta,phi) = (" << beta << ", " << bphi << ")" << endl
	   << "\t\tdeta = " << deta << " local position at 0 " << lEdge1 
	//<< "    \tpith (top,center,bottom) = " << topPitch << " " << pitch << " " << bottomPitch << ", dEta = " << deta 
	//<< ", dPhi = " << dphi 
	   << endl;
  }
  ofos << dashedLine_ << " end" << endl;
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ME0GeometryAnalyzer);

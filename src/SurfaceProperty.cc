#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"

using namespace CLHEP;

static const unsigned flatentries = 2;
static const double minenergy     = 1.0*eV;
static const double maxenergy     = 8.0*eV;


G4OpticalSurface*
MakeS_TyvekCrystal()
{
  const unsigned num       = 2;
  double Ephoton[num]      = {1.5*eV, 8.0*eV};
  double Reflectivity[num] = {0.979, 0.979};

  //////////////////////////////////
  // Realistic Crystal-Tyvek surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface( "TyvekOpSurface" );
  surface->SetType( dielectric_LUT );
  surface->SetModel( LUT );
  surface->SetFinish( polishedtyvekair );

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();

  table->AddProperty( "REFLECTIVITY", Ephoton, Reflectivity, num );

  surface->SetMaterialPropertiesTable( table );

  return surface;
}

G4OpticalSurface*
MakeS_ESR()
{
  const unsigned num       = 2;
  double Ephoton[num]      = {1.5*eV, 8.0*eV};
  double Reflectivity[num] = {0.985, 0.985};
  // source: https://www.osti.gov/servlets/purl/1184400
  //////////////////////////////////
  // ESR surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface( "ESROpSurface" );
  surface->SetType( dielectric_LUT );
  surface->SetModel( LUT );
  surface->SetFinish( polishedvm2000air );


  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();

  table->AddProperty( "REFLECTIVITY", Ephoton, Reflectivity, num );

  surface->SetMaterialPropertiesTable( table );

  return surface;
}

G4OpticalSurface*
MakeS_IdealTyvekCrystal()
{
  //////////////////////////////////
  // Ideal Crystal-Tyvek surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface( "IdealTyvekOpSurface" );
  surface->SetType( dielectric_LUT );
  surface->SetModel( LUT );
  surface->SetFinish( polishedtyvekair );

  return surface;
}
/*
G4OpticalSurface*
MakeS_UnifiedTyvek()
{
  const unsigned nentries        = 2;
  double phoE[nentries]          = {2.0*eV, 6.0*eV};
  double specularlobe[nentries]  = {1.0, 1.0};
  double specularspike[nentries] = {0., 0.};
  double backscatter[nentries]   = {0., 0.};
  double rindex[nentries]        = {1.0, 1.0};
  double reflectivity[nentries]  = {0.979, 0.979};
  double efficiency[nentries]    = {0.0, 0.0};

  //////////////////////////////////
  // Unified Tyvek surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface( "UnifiedTyvekOpSurface" );
  surface->SetType( dielectric_dielectric );
  surface->SetModel( unified );
  surface->SetFinish( groundbackpainted );
  surface->SetSigmaAlpha( 1.3*degree );

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();

  table->AddProperty( "RINDEX",                phoE, rindex,        nentries );
  table->AddProperty( "SPECULARLOBECONSTANT",  phoE, specularlobe,  nentries );
  table->AddProperty( "SPECULARSPIKECONSTANT", phoE, specularspike, nentries );
  table->AddProperty( "BACKSCATTERCONSTANT",   phoE, backscatter,   nentries );
  table->AddProperty( "REFLECTIVITY",          phoE, reflectivity,  nentries );
  table->AddProperty( "EFFICIENCY",            phoE, efficiency,    nentries );

  surface->SetMaterialPropertiesTable( table );

  return surface;
}

G4OpticalSurface*
MakeS_IdealTyvek()
{
  const unsigned nentries        = 2;
  double phoE[nentries]          = {2.0*eV, 6.0*eV};
  double specularlobe[nentries]  = {1.0, 1.0};
  double specularspike[nentries] = {0., 0.};
  double backscatter[nentries]   = {0., 0.};
  double rindex[nentries]        = {1.0, 1.0};
  double reflectivity[nentries]  = {1.0, 1.0};
  double efficiency[nentries]    = {0.0, 0.0};

  //////////////////////////////////
  // Unified Ideal Tyvek surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface(
    "UnifiedIdealTyvekOpSurface" );
  surface->SetType( dielectric_dielectric );
  surface->SetModel( unified );
  surface->SetFinish( groundbackpainted );
  surface->SetSigmaAlpha( 1.3*degree );

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();

  table->AddProperty( "RINDEX",                phoE, rindex,        nentries );
  table->AddProperty( "SPECULARLOBECONSTANT",  phoE, specularlobe,  nentries );
  table->AddProperty( "SPECULARSPIKECONSTANT", phoE, specularspike, nentries );
  table->AddProperty( "BACKSCATTERCONSTANT",   phoE, backscatter,   nentries );
  table->AddProperty( "REFLECTIVITY",          phoE, reflectivity,  nentries );
  table->AddProperty( "EFFICIENCY",            phoE, efficiency,    nentries );

  surface->SetMaterialPropertiesTable( table );

  return surface;
}
*/
G4OpticalSurface*
MakeS_Polished()
{
  static const unsigned nentries       = 2;
  static double phoE[nentries]         = {1.0*eV, 6.0*eV};
  static double specularlobe[nentries] = {1.0, 1.0};

  //////////////////////////////////
  // Realistic polished surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface( "PolishedOpSurface" );
  surface->SetType( dielectric_dielectric );
  surface->SetModel( unified );
  surface->SetFinish( ground );
  // necessary even for polished surfaces to enable UNIFIED code
  surface->SetSigmaAlpha( 1.3 * degree );// Janecek2010 (1.3 * degree)

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();

  table->AddProperty(
    "SPECULARLOBECONSTANT", phoE, specularlobe, nentries );
  surface->SetMaterialPropertiesTable( table );

  return surface;
}

G4OpticalSurface*
MakeS_IdealPolished()
{
  //////////////////////////////////
  // Ideal polished surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface( "IdealOpSurface" );
  surface->SetType( dielectric_dielectric );
  surface->SetModel( glisur );
  surface->SetFinish( polished );

  return surface;
}

G4OpticalSurface*
MakeS_Mirror()
{
  const unsigned nentries       = 2;
  double phoE[nentries]         = {1.0*eV, 6.0*eV};
  double reflectivity[nentries] = {0.9, 0.9};
  //////////////////////////////////
  // Mirror surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface( "MirrorOpSurface" );
  surface->SetType( dielectric_metal );
  surface->SetFinish( polished );
  surface->SetModel( unified );

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();
  table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  surface->SetMaterialPropertiesTable( table );

  return surface;
}

G4OpticalSurface*
MakeS_IdealMirror()
{
  const unsigned nentries       = 2;
  double phoE[nentries]         = {1.0*eV, 6.0*eV};
  double reflectivity[nentries] = {1.0, 1.0};
  //////////////////////////////////
  // Ideal mirror surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface( "MirrorOpSurface" );
  surface->SetType( dielectric_metal );
  surface->SetFinish( polished );
  surface->SetModel( unified );

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();
  table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  surface->SetMaterialPropertiesTable( table );

  return surface;
}

G4OpticalSurface*
MakeS_IdealWhiteSurface()
{
  //////////////////////////////////
  // Ideal mirror surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface( "WhiteOpSurface" );
  surface->SetType( dielectric_metal );
  surface->SetFinish( ground );
  surface->SetModel( unified );

  double phoE[flatentries]         = {1*eV, 6*eV};
  double reflectivity[flatentries] = {0.0, 0.0};

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();
  table->AddProperty( "REFLECTIVITY", phoE, reflectivity, flatentries );
  surface->SetMaterialPropertiesTable( table );

  return surface;
}

G4OpticalSurface*
MakeS_Absorbing()
{
  const unsigned nentries       = 2;
  double phoE[nentries]         = {1.0*eV, 6.0*eV};
  double reflectivity[nentries] = {0.0, 0.0};
  //////////////////////////////////
  // Absorbing surface
  //////////////////////////////////
  G4OpticalSurface* surface = new G4OpticalSurface( "AbsorbingOpSurface" );
  surface->SetType( dielectric_dielectric );
  surface->SetFinish( groundfrontpainted );
  surface->SetModel( unified );

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();
  table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  surface->SetMaterialPropertiesTable( table );

  return surface;
}

G4OpticalSurface*
MakeS_wrap()
{
  //////////////////////////////////
  // wrap surface
  //////////////////////////////////
    G4OpticalSurface* scintWrap = new G4OpticalSurface("ScintWrap");
    scintWrap->SetType(dielectric_metal);
    scintWrap->SetFinish(ground);
    scintWrap->SetModel(unified);
    double pp[2]={0.1*eV,6*eV};
    double reflectivity[2] = {0.05, 0.05};
    G4MaterialPropertiesTable* scintWrapProperty
      = new G4MaterialPropertiesTable();
    scintWrapProperty->AddProperty("REFLECTIVITY",pp,reflectivity,2);
    scintWrap->SetMaterialPropertiesTable(scintWrapProperty);
  return scintWrap;
}

G4OpticalSurface* MakeS_PMT(){
  //surface of detector
  const unsigned nentries = 21;
  double phoE[nentries]   = { 
4.591666802	*eV,
4.249963562	*eV,
3.970289068	*eV,
3.545647077	*eV,
3.157143617	*eV,
2.8391158	*eV,
2.630406004	*eV,
2.436700346	*eV,
2.337427879	*eV,
2.246008967	*eV,
2.193827957	*eV,
2.102995234	*eV,
2.041160762	*eV,
1.99177784	*eV,
1.910653049	*eV,
1.872420367	*eV,
1.830538262	*eV,
1.812740321	*eV,
1.78790683	*eV,
1.768575074	*eV,
1.761389217     *eV	
};
  double efficiency[nentries] = {
13.033	,
20.201	,
23.97	,
25.352	,
25.317	,
22.981	,
18.968	,
15.068	,
9.168	,
6.3757	,
5.1664	,
3.3926	,
2.2282	,
1.3305	,
0.40721	,
0.18973	,
0.073042	,
0.043624	,
0.021119	,
0.012856	,
0.010031	
};
for(int dd=0;dd < 21;dd++) efficiency[dd]/=100.;

  const unsigned ref_ent = 2;
  double phoE2[ref_ent]   = {1*eV, 6*eV};
  double refraction[ref_ent] = {1.473, 1.473};
  double QE[ref_ent] = {1.0, 1.0};

  G4OpticalSurface* surface = new G4OpticalSurface( "PMT_Surface" );
  surface->SetType( dielectric_metal );
  surface->SetModel( unified );
//use 100% QE
  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();
  table->AddProperty( "EFFICIENCY",   phoE2,  QE,   ref_ent );

//https://www.phys.hawaii.edu/~varner/OpticaConstants_K2CsSb_RbCsSb_2005.pdf
  double reflectivity[ref_ent] = {0.16, 0.16};


  table->AddProperty( "REFLECTIVITY", phoE2, reflectivity, ref_ent );
  surface->SetMaterialPropertiesTable( table );



  return surface;

}
/*
G4OpticalSurface* MakeS_SiPM()
{
  // From Meeting on HB SiPM selection (James Hirschauer)
*/
/*
  const unsigned nentries = 49;
  double phoE[nentries]   = {
    3.269640158*eV, 3.191037408*eV, 3.128716912*eV, 3.102662139*eV,
    3.069524008*eV, 3.034012326*eV, 3.001068157*eV, 2.959583079*eV,
    2.901064107*eV, 2.838801776*eV, 2.779808831*eV, 2.721566033*eV,
    2.667612494*eV, 2.615095968*eV, 2.563388259*eV, 2.515844086*eV,
    2.469260675*eV, 2.423979937*eV, 2.379744327*eV, 2.337228676*eV,
    2.296532193*eV, 2.257220927*eV, 2.219239034*eV, 2.182957454*eV,
    2.147148186*eV, 2.112470008*eV, 2.079619445*eV, 2.047338607*eV,
    2.015906691*eV, 1.987528043*eV, 1.959370034*eV, 1.930211969*eV,
    1.902680667*eV, 1.875153451*eV, 1.84878577*eV,  1.823292126*eV,
    1.799066012*eV, 1.775177229*eV, 1.751266161*eV, 1.728338166*eV,
    1.705491878*eV, 1.683828975*eV, 1.662853866*eV, 1.643031375*eV,
    1.623272564*eV, 1.603516479*eV, 1.584163112*eV, 1.565103929*eV,
    1.548241716*eV};

  double efficiency[nentries] = {
    0.00080151,  0.003873974,  0.009957452,   0.025682901,
    0.070237974, 0.150458664,  0.218876073,   0.239857707,
    0.250435989, 0.259508637,  0.264540183,   0.266889254,
    0.270669968, 0.273577562,  0.272402437,   0.270462795,
    0.268157907, 0.264243881,  0.25831461,    0.252029247,
    0.245123015, 0.237845409,  0.230874024,   0.224981219,
    0.220313658, 0.21616835,   0.211206974,   0.205552053,
    0.199704762, 0.193696557,  0.187884494,   0.182241421,
    0.176095301, 0.168483378,  0.16291645,    0.157493577,
    0.152319279, 0.146859841,  0.141516867,   0.135704067,
    0.129615124, 0.123742916,  0.11736541,    0.11184209,
    0.105217764, 0.097865669,  0.091000102,   0.084681686,
    0.080451249};
*/
/*
  const unsigned nentries = 51;
  double phoE[nentries]   = {
4.456977245*eV,
4.405663883*eV,
4.371335649*eV,
4.337993527*eV,
4.305156186*eV,
4.272812248*eV,
4.120034327*eV,
4.047934735*eV,
3.885068561*eV,
3.795976762*eV,
3.710990512*eV,
3.641025285*eV,
3.562661791*eV,
3.416389546*eV,
3.254178294*eV,
3.200252775*eV,
3.130913965*eV,
3.072870849*eV,
3.016867241*eV,
2.948003733*eV,
2.868012792*eV,
2.772393126*eV,
2.622893865*eV,
2.488842802*eV,
2.467789116*eV,
2.382432947*eV,
2.325721122*eV,
2.276066914*eV,
2.224250888*eV,
2.174741594*eV,
2.115879533*eV,
2.067333517*eV,
2.024496146*eV,
1.98339801*eV,
1.943935293*eV,
1.906041585*eV,
1.869568783*eV,
1.834465614*eV,
1.80065635*eV,
1.76540215*eV,
1.731477711*eV,
1.696391875*eV,
1.66504429*eV,
1.628029216*eV,
1.594775069*eV,
1.562852229*eV,
1.546342471*eV,
1.512389673*eV,
1.477990547*eV,
1.443354983*eV,
1.40859115*eV};
  double efficiency[nentries] = {
1.1045,
2.4287,
4.4945,
5.8187,
7.1429,
8.4672,
11.9633,
13.2347,
15.5127,
16.7841,
18.0556,
19.327,
20.5984,
21.6053,
22.5064,
23.7778,
25.7379,
27.0093,
28.2808,
29.4993,
30.7179,
31.7778,
32.0435,
31.7796,
31.6738,
30.7211,
29.5563,
28.3914,
27.1736,
26.0618,
24.6852,
23.5734,
22.4086,
21.2438,
20.079,
18.9142,
17.7494,
16.5846,
15.4198,
14.255,
13.1432,
12.0315,
10.9197,
9.755,
8.6962,
7.6375,
7.1611,
6.1554,
5.3086,
4.4618,
3.6681};
for(int dd=0;dd < 51;dd++) efficiency[dd]/=100.;
  // SiPM reflectivity paper: Using a flat 0.12 for now
  // https://reader.elsevier.com/reader/sd/pii/S016890021831297X?token=FADE142C43F500DEE5289D65B9D6573497F4AD220F3462EBA9239256CFD8E11FF98FF84E302016EAB89FA1B62B1EEDCE
//  double phoE2[flatentries]        = {minenergy, maxenergy };
//  double reflectivity[flatentries] = {0.12, 0.12};

  const unsigned ref_ent = 32;
  double phoE2[ref_ent]   = {
4.928613174	*eV,
4.760932071	*eV,
4.57675131	*eV,
4.416020551	*eV,
4.27782469	*eV,
4.132944198	*eV,
3.988297134	*eV,
3.870877084	*eV,
3.751866883	*eV,
3.647774074	*eV,
3.538260695	*eV,
3.441997529	*eV,
3.347666946	*eV,
3.264460058	*eV,
3.176475533	*eV,
3.098752668	*eV,
3.019365195	*eV,
2.751902006	*eV,
2.432636667	*eV,
2.254544997	*eV,
2.136148464	*eV,
2.032627719	*eV,
1.907273067	*eV,
1.796976535	*eV,
1.722050515	*eV,
1.570692624	*eV,
1.511965476	*eV,
1.459049532	*eV,
1.377862407	*eV,
1.305495288	*eV,
1.265726027	*eV,
1.240350474	*eV
};
double reflectivity[ref_ent] = {
0.67091	,
0.70309	,
0.73346	,
0.72132	,
0.68334	,
0.62105	,
0.59036	,
0.57232	,
0.5645	,
0.56008	,
0.56441	,
0.58236	,
0.5828	,
0.54482	,
0.51025	,
0.48443	,
0.46541	,
0.41855	,
0.38169	,
0.3657	,
0.35705	,
0.35096	,
0.34397	,
0.33772	,
0.33321	,
0.32645	,
0.32498	,
0.32241	,
0.31988	,
0.31625	,
0.31417	,
0.31323	
};


  G4OpticalSurface* surface = new G4OpticalSurface( "SiPM_Surface" );
  surface->SetType( dielectric_metal );
  surface->SetModel( unified );

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();

  table->AddProperty( "EFFICIENCY",   phoE,  efficiency,   nentries );
  table->AddProperty( "REFLECTIVITY", phoE2, reflectivity, ref_ent );
  surface->SetMaterialPropertiesTable( table );
  return surface;
}
*/
G4OpticalSurface*
MakeS_PCBSurface()
{
  // PCB is a flat gray surface for now
  double phoE2[flatentries]        = {minenergy, maxenergy };
  double reflectivity[flatentries] = {0.5, 0.5};

  G4OpticalSurface* surface = new G4OpticalSurface( "PCB_Surface" );
  surface->SetType( dielectric_metal );
  surface->SetFinish( ground );
  surface->SetModel( unified );

  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable();

  table->AddProperty( "REFLECTIVITY", phoE2, reflectivity, flatentries );
  surface->SetMaterialPropertiesTable( table );
  return surface;
}


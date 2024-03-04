######################################################################
#
#  custom steering file for ODD simulation 
#  It applies a virtual scoring of energy around the incident particle
#
#
######################################################################
from DDSim.DD4hepSimulation import DD4hepSimulation
import DDG4
from g4units import m, mm, GeV, MeV, rad, keV,TeV
import os
from math import pi

SIM = DD4hepSimulation()

## The compact XML file
if "FCCcore" in os.environ:
   SIM.compactFile = os.environ["FCCcore"]+"/OpenDataDetector/xml/OpenDataDetector_noTracker.xml"
   print(f"Detector: {os.environ['FCCcore']}/OpenDataDetector/xml/OpenDataDetector_noTracker.xml")
## Lorentz boost for the crossing angle, in radian!
SIM.enableDetailedShowerMode = False
SIM.enableG4GPS = True
SIM.enableG4Gun = False
SIM.enableGun = False
## InputFiles for simulation .stdhep, .slcio, .HEPEvt, .hepevt, .hepmc, .pairs files are supported
#SIM.inputFiles = [""]
## Macro file to execute for runType 'run' or 'vis'
if "FCCcore" in os.environ:
   SIM.macroFile = os.environ["FCCcore"]+"/ddfastsim/options/ODD/ECalBarrelEndcap_singleEmMesh.mac"
   print(f"Detector: {os.environ['FCCcore']}/ddfastsim/options/ODD/ECalBarrelEndcap_singleEmMesh.mac")
## number of events to simulate, used in batch mode
SIM.numberOfEvents = 100
## Outputfile from the simulation,only lcio output is supported
#SIM.outputFile = "dummyOutput_edm4hep.root" ##"dummyOutput.slcio"
SIM.outputFile = "singleShower_Par04_lowgran_1GeVto1TeVpower_100showers_gamma_thetaUserCone_phiFull_edm4hep.root"

## Physics list to use in simulation
SIM.physicsList = None
## Verbosity use integers from 1(most) to 7(least) verbose
## or strings: VERBOSE, DEBUG, INFO, WARNING, ERROR, FATAL, ALWAYS
SIM.printLevel = "INFO"
## The type of action to do in this invocation
## batch: just simulate some events, needs numberOfEvents, and input file or gun
## vis: enable visualisation, run the macroFile if it is set
## run: run the macroFile and exit
## shell: enable interactive session
SIM.runType = "run" # "batch"
## Skip first N events when reading a file
SIM.skipNEvents = 0
## Steering file to change default behaviour
SIM.steeringFile = None
## FourVector of translation for the Smearing of the Vertex position: x y z t
SIM.vertexOffset = [0.0, 0.0, 0.0, 0.0]
## FourVector of the Sigma for the Smearing of the Vertex position: x y z t
SIM.vertexSigma = [0.0, 0.0, 0.0, 0.0]


################################################################################
## Action holding sensitive detector actions
##   The default tracker and calorimeter actions can be set with
##
##   >>> SIM = DD4hepSimulation()
##   >>> SIM.action.tracker = ('Geant4TrackerWeightedAction', {'HitPositionCombination': 2, 'CollectSingleDeposits': False})
##   >>> SIM.action.calo    = "Geant4CalorimeterAction"
##
##   for specific subdetectors specific sensitive detectors can be set based on pattern matching
##
##   >>> SIM = DD4hepSimulation()
##   >>> SIM.action.mapActions['tpc'] = "TPCSDAction"
##
##   and additional parameters for the sensitive detectors can be set when the map is given a tuple
##
##   >>> SIM = DD4hepSimulation()
##   >>> SIM.action.mapActions['ecal'] =( "CaloPreShowerSDAction", {"FirstLayerNumber": 1} )
##
##
################################################################################

from math import pi
meshCellNumRho = 9
meshCellNumPhi = 16
meshCellNumZ = 45
meshCellSizeRho = 4.65 * mm
meshCellSizePhi = 2 * pi / meshCellNumPhi
meshCellSizeZ = 5.05 * mm


##  set the default calorimeter action
SIM.action.calo = ('Geant4ShowerMeshCalorimeterAction', {'sizeOfZCells': meshCellSizeZ,
'sizeOfRhoCells': meshCellSizeRho, 'sizeOfPhiCells': meshCellSizePhi,
'numOfZCells': meshCellNumZ, 'numOfRhoCells': meshCellNumRho, })

##  create a map of patterns and actions to be applied to sensitive detectors
##         example: SIM.action.mapActions['tpc'] = "TPCSDAction"
SIM.action.mapActions = {}

##  set the default tracker action 
SIM.action.tracker =  ('Geant4TrackerWeightedAction', {'HitPositionCombination': 2, 'CollectSingleDeposits': False})

# Configure Run actions
kernel = DDG4.Kernel()
run1 = DDG4.RunAction(kernel,'ScoreMeshRunAction/runaction')
run1.sizeOfZCells= meshCellSizeZ
run1.sizeOfRhoCells= meshCellSizeRho
run1.sizeOfPhiCells= meshCellSizePhi
run1.numOfZCells= meshCellNumZ
run1.numOfRhoCells= meshCellNumRho
kernel.registerGlobalAction(run1)
kernel.runAction().add(run1)
event1 = DDG4.EventAction(kernel,'ScoreMeshEventAction/eventaction')
event1.sizeOfZCells= meshCellSizeZ
event1.sizeOfRhoCells= meshCellSizeRho
event1.sizeOfPhiCells= meshCellSizePhi
event1.numOfZCells= meshCellNumZ
event1.numOfRhoCells= meshCellNumRho
kernel.registerGlobalAction(event1)
kernel.eventAction().add(event1)


################################################################################
## Configuration for the magnetic field (stepper)
################################################################################
## --- used in v01-19-05 :
#SIM.field.delta_chord = 1e-05
#SIM.field.delta_intersection = 1e-05
#SIM.field.delta_one_step = .5e-03*mm
#SIM.field.eps_max = 1e-04
#SIM.field.eps_min = 1e-05
#SIM.field.equation = "Mag_UsualEqRhs"
#SIM.field.largest_step = 10.*m
#SIM.field.min_chord_step = 1.e-2*mm
#SIM.field.stepper = "HelixSimpleRunge"

################################################################################
## Configuration for sensitive detector filters
##
##   Set the default filter for tracker or caliromter
##   >>> SIM.filter.tracker = "edep1kev"
##   >>> SIM.filter.calo = ""
##
##   Assign a filter to a sensitive detector via pattern matching
##   >>> SIM.filter.mapDetFilter['FTD'] = "edep1kev"
##
##   Or more than one filter:
##   >>> SIM.filter.mapDetFilter['FTD'] = ["edep1kev", "geantino"]
##
##   Don't use the default filter or anything else:
##   >>> SIM.filter.mapDetFilter['TPC'] = None ## or "" or []
##
##   Create a custom filter. The dictionary is used to instantiate the filter later on
##   >>> SIM.filter.filters['edep3kev'] = dict(name="EnergyDepositMinimumCut/3keV", parameter={"Cut": 3.0*keV} )
##
##
################################################################################

##  default filter for calorimeter sensitive detectors; this is applied if no other filter is used for a calorimeter
SIM.filter.calo = "edep0"

##  list of filter objects: map between name and parameter dictionary
SIM.filter.filters = {'edep0': {'parameter': {'Cut': 0.0}, 'name': 'EnergyDepositMinimumCut/Cut0'}, 'geantino': {'parameter': {}, 'name': 'GeantinoRejectFilter/GeantinoRejector'}, 'edep1kev': {'parameter': {'Cut': 0.001}, 'name': 'EnergyDepositMinimumCut'}}

##  a map between patterns and filter objects, using patterns to attach filters to sensitive detector
SIM.filter.mapDetFilter = {}

# ##  default filter for tracking sensitive detectors; this is applied if no other filter is used for a tracker
# SIM.filter.tracker = "edep0"


################################################################################
## Configuration for the GuineaPig InputFiles
################################################################################

## Set the number of pair particles to simulate per event.
##     Only used if inputFile ends with ".pairs"
##     If "-1" all particles will be simulated in a single event
##
SIM.guineapig.particlesPerEvent = "-1"


################################################################################
## Configuration for the DDG4 ParticleGun
################################################################################


##  direction of the particle gun, 3 vector
SIM.gun.direction = (0, 1, 0)

## choose the distribution of the random direction for theta
##
##     Options for random distributions:
##
##     'uniform' is the default distribution, flat in theta
##     'cos(theta)' is flat in cos(theta)
##     'eta', or 'pseudorapidity' is flat in pseudorapity
##     'ffbar' is distributed according to 1+cos^2(theta)
##
##     Setting a distribution will set isotrop = True
##
SIM.gun.distribution = "uniform"

## Total energy (including mass) for the particle gun.
##
## If not None, it will overwrite the setting of momentumMin and momentumMax
SIM.gun.energy = 0

##  isotropic distribution for the particle gun
##
##     use the options phiMin, phiMax, thetaMin, and thetaMax to limit the range of randomly distributed directions
##     if one of these options is not None the random distribution will be set to True and cannot be turned off!
##
#SIM.gun.isotrop = True

## Maximal momentum when using distribution (default = 0.0)
SIM.gun.momentumMax = 10*GeV

## Minimal momentum when using distribution (default = 0.0)
SIM.gun.momentumMin = 1*GeV
SIM.gun.multiplicity = 1
SIM.gun.particle = "gamma"

## Maximal azimuthal angle for random distribution
SIM.gun.phiMax = 0

## Minimal azimuthal angle for random distribution
SIM.gun.phiMin = 0

##  position of the particle gun, 3 vector
SIM.gun.position = (0.0, 0.0, 0.0)

## Maximal polar angle for random distribution
#SIM.gun.thetaMax = 0

## Minimal polar angle for random distribution
#SIM.gun.thetaMin = 0

## Maximal pseudorapidity for random distibution (overrides thetaMin)
SIM.gun.etaMax = 0

## Minimal pseudorapidity for random distibution (overrides thetaMax)
SIM.gun.etaMin = 0

################################################################################
## Configuration for the output levels of DDG4 components
################################################################################

## Output level for input sources
SIM.output.inputStage = 3

## Output level for Geant4 kernel
SIM.output.kernel = 3

## Output level for ParticleHandler
SIM.output.part = 3

## Output level for Random Number Generator setup
SIM.output.random = 6


################################################################################
## Configuration for the Particle Handler/ MCTruth treatment
################################################################################

## Enable lots of printout on simulated hits and MC-truth information
SIM.part.enableDetailedHitsAndParticleInfo = False

##  Keep all created particles
SIM.part.keepAllParticles = False

## Minimal distance between particle vertex and endpoint of parent after
##     which the vertexIsNotEndpointOfParent flag is set
##
SIM.part.minDistToParentVertex = 2.2e-14

## MinimalKineticEnergy to store particles created in the tracking region
# This setting vastly reduces the file size as only the incident particle is stored
SIM.part.minimalKineticEnergy = 1*TeV

##  Printout at End of Tracking
SIM.part.printEndTracking = False

##  Printout at Start of Tracking
SIM.part.printStartTracking = False

## List of processes to save, on command line give as whitespace separated string in quotation marks
SIM.part.saveProcesses = []

## Optionally enable an extended Particle Handler
SIM.part.userParticleHandler = ""

################################################################################
## Configuration for the PhysicsList
################################################################################
# this needs to be set to False if any standard physics list is used:
SIM.physics.decays = False
SIM.physics.list = "FTFP_BERT"

##  location of particle.tbl file containing extra particles and their lifetime information
##
SIM.physics.pdgfile = os.path.join( os.environ.get("DD4HEP"),  "DDG4/examples/particle.tbl")

##  The global geant4 rangecut for secondary production
##
##     Default is 0.7 mm as is the case in geant4 10
##
##     To disable this plugin and be absolutely sure to use the Geant4 default range cut use "None"
##
##     Set printlevel to DEBUG to see a printout of all range cuts,
##     but this only works if range cut is not "None"
##
SIM.physics.rangecut =  0.1*mm


################################################################################
## Properties for the random number generator
################################################################################

## If True, calculate random seed for each event based on eventID and runID
## allows reproducibility even when SkippingEvents
SIM.random.enableEventSeed = False
SIM.random.file = None
SIM.random.luxury = 1
SIM.random.replace_gRandom = True
SIM.random.seed = None
SIM.random.type = None

################################################################################
## Configure fast simulation
################################################################################

def fastsimSettings(kernel):
   from g4units import GeV, MeV  # DO NOT REMOVE OR MOVE!!!!! (EXCLAMATION MARK)
   from DDG4 import DetectorConstruction, Geant4, PhysicsList

   geant4 = Geant4(kernel)
   seq = geant4.detectorConstruction()

   # Create a model for fast simulation
   model = DetectorConstruction(kernel, "DefineMeshModel" )
   # Mandatory model parameters
   model.RegionName = 'ECalBarrelRegion'
   model.Enable = True
   allParticles = [ "B+", "B-", "B0", "Bc+", "Bc-", "Bs0", "D+", "D-", "D0", "Ds+", "Ds-", "GenericIon", "He3", "J/psi", "N(1440)+", "N(1440)0", "N(1520)+", "N(1520)0", "N(1535)+", "N(1535)0", "N(1650)+", "N(1650)0", "N(1675)+", "N(1675)0", "N(1680)+", "N(1680)0", "N(1700)+", "N(1700)0", "N(1710)+", "N(1710)0", "N(1720)+", "N(1720)0", "N(1900)+", "N(1900)0", "N(1990)+", "N(1990)0", "N(2090)+", "N(2090)0", "N(2190)+", "N(2190)0", "N(2220)+", "N(2220)0", "N(2250)+", "N(2250)0", "Upsilon", "a0(1450)+", "a0(1450)-", "a0(1450)0", "a0(980)+", "a0(980)-", "a0(980)0", "a1(1260)+", "a1(1260)-", "a1(1260)0", "a2(1320)+", "a2(1320)-", "a2(1320)0", "alpha", "anti_B0", "anti_Bs0", "anti_D0", "anti_He3", "anti_N(1440)+", "anti_N(1440)0", "anti_N(1520)+", "anti_N(1520)0", "anti_N(1535)+", "anti_N(1535)0", "anti_N(1650)+", "anti_N(1650)0", "anti_N(1675)+", "anti_N(1675)0", "anti_N(1680)+", "anti_N(1680)0", "anti_N(1700)+", "anti_N(1700)0", "anti_N(1710)+", "anti_N(1710)0", "anti_N(1720)+", "anti_N(1720)0", "anti_N(1900)+", "anti_N(1900)0", "anti_N(1990)+", "anti_N(1990)0", "anti_N(2090)+", "anti_N(2090)0", "anti_N(2190)+", "anti_N(2190)0", "anti_N(2220)+", "anti_N(2220)0", "anti_N(2250)+", "anti_N(2250)0", "anti_alpha", "anti_b_quark", "anti_bb1_diquark", "anti_bc0_diquark", "anti_bc1_diquark", "anti_bd0_diquark", "anti_bd1_diquark", "anti_bs0_diquark", "anti_bs1_diquark", "anti_bu0_diquark", "anti_bu1_diquark", "anti_c_quark", "anti_cc1_diquark", "anti_cd0_diquark", "anti_cd1_diquark", "anti_cs0_diquark", "anti_cs1_diquark", "anti_cu0_diquark", "anti_cu1_diquark", "anti_d_quark", "anti_dd1_diquark", "anti_delta(1600)+", "anti_delta(1600)++", "anti_delta(1600)-", "anti_delta(1600)0", "anti_delta(1620)+", "anti_delta(1620)++", "anti_delta(1620)-", "anti_delta(1620)0", "anti_delta(1700)+", "anti_delta(1700)++", "anti_delta(1700)-", "anti_delta(1700)0", "anti_delta(1900)+", "anti_delta(1900)++", "anti_delta(1900)-", "anti_delta(1900)0", "anti_delta(1905)+", "anti_delta(1905)++", "anti_delta(1905)-", "anti_delta(1905)0", "anti_delta(1910)+", "anti_delta(1910)++", "anti_delta(1910)-", "anti_delta(1910)0", "anti_delta(1920)+", "anti_delta(1920)++", "anti_delta(1920)-", "anti_delta(1920)0", "anti_delta(1930)+", "anti_delta(1930)++", "anti_delta(1930)-", "anti_delta(1930)0", "anti_delta(1950)+", "anti_delta(1950)++", "anti_delta(1950)-", "anti_delta(1950)0", "anti_delta+", "anti_delta++", "anti_delta-", "anti_delta0", "anti_deuteron", "anti_doublehyperH4", "anti_doublehyperdoubleneutron", "anti_hyperH4", "anti_hyperHe5", "anti_hyperalpha", "anti_hypertriton", "anti_k(1460)0", "anti_k0_star(1430)0", "anti_k1(1270)0", "anti_k1(1400)0", "anti_k2(1770)0", "anti_k2_star(1430)0", "anti_k2_star(1980)0", "anti_k3_star(1780)0", "anti_k_star(1410)0", "anti_k_star(1680)0", "anti_k_star0", "anti_kaon0", "anti_lambda", "anti_lambda(1405)", "anti_lambda(1520)", "anti_lambda(1600)", "anti_lambda(1670)", "anti_lambda(1690)", "anti_lambda(1800)", "anti_lambda(1810)", "anti_lambda(1820)", "anti_lambda(1830)", "anti_lambda(1890)", "anti_lambda(2100)", "anti_lambda(2110)", "anti_lambda_b", "anti_lambda_c+", "anti_neutron", "anti_nu_e", "anti_nu_mu", "anti_nu_tau", "anti_omega-", "anti_omega_b-", "anti_omega_c0", "anti_proton", "anti_s_quark", "anti_sd0_diquark", "anti_sd1_diquark", "anti_sigma(1385)+", "anti_sigma(1385)-", "anti_sigma(1385)0", "anti_sigma(1660)+", "anti_sigma(1660)-", "anti_sigma(1660)0", "anti_sigma(1670)+", "anti_sigma(1670)-", "anti_sigma(1670)0", "anti_sigma(1750)+", "anti_sigma(1750)-", "anti_sigma(1750)0", "anti_sigma(1775)+", "anti_sigma(1775)-", "anti_sigma(1775)0", "anti_sigma(1915)+", "anti_sigma(1915)-", "anti_sigma(1915)0", "anti_sigma(1940)+", "anti_sigma(1940)-", "anti_sigma(1940)0", "anti_sigma(2030)+", "anti_sigma(2030)-", "anti_sigma(2030)0", "anti_sigma+", "anti_sigma-", "anti_sigma0", "anti_sigma_b+", "anti_sigma_b-", "anti_sigma_b0", "anti_sigma_c+", "anti_sigma_c++", "anti_sigma_c0", "anti_ss1_diquark", "anti_su0_diquark", "anti_su1_diquark", "anti_t_quark", "anti_triton", "anti_u_quark", "anti_ud0_diquark", "anti_ud1_diquark", "anti_uu1_diquark", "anti_xi(1530)-", "anti_xi(1530)0", "anti_xi(1690)-", "anti_xi(1690)0", "anti_xi(1820)-", "anti_xi(1820)0", "anti_xi(1950)-", "anti_xi(1950)0", "anti_xi(2030)-", "anti_xi(2030)0", "anti_xi-", "anti_xi0", "anti_xi_b-", "anti_xi_b0", "anti_xi_c+", "anti_xi_c0", "b1(1235)+", "b1(1235)-", "b1(1235)0", "b_quark", "bb1_diquark", "bc0_diquark", "bc1_diquark", "bd0_diquark", "bd1_diquark", "bs0_diquark", "bs1_diquark", "bu0_diquark", "bu1_diquark", "c_quark", "cc1_diquark", "cd0_diquark", "cd1_diquark", "chargedgeantino", "cs0_diquark", "cs1_diquark", "cu0_diquark", "cu1_diquark", "d_quark", "dd1_diquark", "delta(1600)+", "delta(1600)++", "delta(1600)-", "delta(1600)0", "delta(1620)+", "delta(1620)++", "delta(1620)-", "delta(1620)0", "delta(1700)+", "delta(1700)++", "delta(1700)-", "delta(1700)0", "delta(1900)+", "delta(1900)++", "delta(1900)-", "delta(1900)0", "delta(1905)+", "delta(1905)++", "delta(1905)-", "delta(1905)0", "delta(1910)+", "delta(1910)++", "delta(1910)-", "delta(1910)0", "delta(1920)+", "delta(1920)++", "delta(1920)-", "delta(1920)0", "delta(1930)+", "delta(1930)++", "delta(1930)-", "delta(1930)0", "delta(1950)+", "delta(1950)++", "delta(1950)-", "delta(1950)0", "delta+", "delta++", "delta-", "delta0", "deuteron", "doublehyperH4", "doublehyperdoubleneutron", "e+", "e-", "eta", "eta(1295)", "eta(1405)", "eta(1475)", "eta2(1645)", "eta2(1870)", "eta_prime", "etac", "f0(1370)", "f0(1500)", "f0(1710)", "f0(500)", "f0(980)", "f1(1285)", "f1(1420)", "f2(1270)", "f2(1810)", "f2(2010)", "f2_prime(1525)", "gamma", "geantino", "gluon", "h1(1170)", "h1(1380)", "hyperH4", "hyperHe5", "hyperalpha", "hypertriton", "k(1460)+", "k(1460)-", "k(1460)0", "k0_star(1430)+", "k0_star(1430)-", "k0_star(1430)0", "k1(1270)+", "k1(1270)-", "k1(1270)0", "k1(1400)+", "k1(1400)-", "k1(1400)0", "k2(1770)+", "k2(1770)-", "k2(1770)0", "k2_star(1430)+", "k2_star(1430)-", "k2_star(1430)0", "k2_star(1980)+", "k2_star(1980)-", "k2_star(1980)0", "k3_star(1780)+", "k3_star(1780)-", "k3_star(1780)0", "k_star(1410)+", "k_star(1410)-", "k_star(1410)0", "k_star(1680)+", "k_star(1680)-", "k_star(1680)0", "k_star+", "k_star-", "k_star0", "kaon+", "kaon-", "kaon0", "kaon0L", "kaon0S", "lambda", "lambda(1405)", "lambda(1520)", "lambda(1600)", "lambda(1670)", "lambda(1690)", "lambda(1800)", "lambda(1810)", "lambda(1820)", "lambda(1830)", "lambda(1890)", "lambda(2100)", "lambda(2110)", "lambda_b", "lambda_c+", "mu+", "mu-", "neutron", "nu_e", "nu_mu", "nu_tau", "omega", "omega(1420)", "omega(1650)", "omega-", "omega3(1670)", "omega_b-", "omega_c0", "opticalphoton", "phi", "phi(1680)", "phi3(1850)", "pi(1300)+", "pi(1300)-", "pi(1300)0", "pi+", "pi-", "pi0", "pi2(1670)+", "pi2(1670)-", "pi2(1670)0", "proton", "rho(1450)+", "rho(1450)-", "rho(1450)0", "rho(1700)+", "rho(1700)-", "rho(1700)0", "rho+", "rho-", "rho0", "rho3(1690)+", "rho3(1690)-", "rho3(1690)0", "s_quark", "sd0_diquark", "sd1_diquark", "sigma(1385)+", "sigma(1385)-", "sigma(1385)0", "sigma(1660)+", "sigma(1660)-", "sigma(1660)0", "sigma(1670)+", "sigma(1670)-", "sigma(1670)0", "sigma(1750)+", "sigma(1750)-", "sigma(1750)0", "sigma(1775)+", "sigma(1775)-", "sigma(1775)0", "sigma(1915)+", "sigma(1915)-", "sigma(1915)0", "sigma(1940)+", "sigma(1940)-", "sigma(1940)0", "sigma(2030)+", "sigma(2030)-", "sigma(2030)0", "sigma+", "sigma-", "sigma0", "sigma_b+", "sigma_b-", "sigma_b0", "sigma_c+", "sigma_c++", "sigma_c0", "ss1_diquark", "su0_diquark", "su1_diquark", "t_quark", "tau+", "tau-", "triton", "u_quark", "ud0_diquark", "ud1_diquark", "uu1_diquark", "xi(1530)-", "xi(1530)0", "xi(1690)-", "xi(1690)0", "xi(1820)-", "xi(1820)0", "xi(1950)-", "xi(1950)0", "xi(2030)-", "xi(2030)0", "xi-", "xi0", "xi_b-", "xi_b0", "xi_c+", "xi_c0"]
   model.ApplicableParticles = allParticles
   model.enableUI()
   seq.adopt(model)
   #-------------------
   model_endcap = DetectorConstruction(kernel, "DefineMeshModel")
   # Mandatory model parameters
   model_endcap.RegionName = 'ECalEndcapRegion'
   model_endcap.Enable = True
   model_endcap.ApplicableParticles = allParticles
   model_endcap.enableUI()
   seq.adopt(model_endcap)
   #-------------------
   modelHcal = DetectorConstruction(kernel, "DefineMeshModel" )
   # Mandatory model parameters
   modelHcal.RegionName = 'HCalBarrelRegion'
   modelHcal.Enable = True
   modelHcal.ApplicableParticles = allParticles
   modelHcal.enableUI()
   seq.adopt(modelHcal)
   #-------------------
   modelHcal_endcap = DetectorConstruction(kernel, "DefineMeshModel")
   # Mandatory model parameters
   modelHcal_endcap.RegionName = 'HCalEndcapRegion'
   modelHcal_endcap.Enable = True
   modelHcal_endcap.ApplicableParticles = allParticles
   modelHcal_endcap.enableUI()
   seq.adopt(modelHcal_endcap)
   #-------------------

   # Now build the physics list:
   phys = kernel.physicsList()
   ph = PhysicsList(kernel, str('Geant4FastPhysics/FastPhysicsList'))
   ph.EnabledParticles = allParticles
   ph.BeVerbose = True
   ph.enableUI()
   phys.adopt(ph)
   phys.dump()

SIM.physics.setupUserPhysics( fastsimSettings)

CXX=g++
OPTIMIZATION=-O0
DEBUG=-g
INCLUDE_DIRS=-I/projects/dcm/ProteinVolumeCalculator
CXXFLAGS=$(OPTIMIZATION) $(DEBUG) $(INCLUDE_DIRS)


	
OBJS =  bmpg_uncc_edu/algorithms/FastGraph.o \
	bmpg_uncc_edu/algorithms/GraphTraits.o \
	bmpg_uncc_edu/algorithms/HashGrid.o \
	bmpg_uncc_edu/chemistry/Atom.o \
	bmpg_uncc_edu/chemistry/Bond.o \
	bmpg_uncc_edu/chemistry/Complex.o \
	bmpg_uncc_edu/chemistry/CovalentBond.o \
	bmpg_uncc_edu/chemistry/CovalentBonder.o \
	bmpg_uncc_edu/chemistry/Element.o \
	bmpg_uncc_edu/chemistry/Molecule.o \
	bmpg_uncc_edu/chemistry/PDBAtom.o \
	bmpg_uncc_edu/chemistry/PDBProteinChain.o \
	bmpg_uncc_edu/chemistry/PDBProteinResidue.o \
	bmpg_uncc_edu/chemistry/detail/PDBChainManager.o \
	bmpg_uncc_edu/chemistry/helper/PDBHelper.o \
	bmpg_uncc_edu/chemistry/library/AminoAcidLibrary.o \
	bmpg_uncc_edu/chemistry/library/BondDistanceTable.o \
	bmpg_uncc_edu/chemistry/library/PeriodicTable.o \
	bmpg_uncc_edu/chemistry/library/loader/AALibLoader.o \
	bmpg_uncc_edu/fast/Defs.o \
	bmpg_uncc_edu/fast/ParameterFile.o \
	bmpg_uncc_edu/fast/input/ResourceBundle.o \
	bmpg_uncc_edu/math/MathConstants.o \
	bmpg_uncc_edu/math/R3Vector.o \
	bmpg_uncc_edu/math/Vector.o \
	bmpg_uncc_edu/util/Exception.o \
	bmpg_uncc_edu/util/FlatFileReader.o \
	bmpg_uncc_edu/util/Properties.o \
	bmpg_uncc_edu/util/SequentialIdGenerator.o \
	bmpg_uncc_edu/util/StringUtil.o \
	bmpg_uncc_edu/util/logger/DefaultLogger.o \
	bmpg_uncc_edu/util/logger/Logger.o \
	bmpg_uncc_edu/util/logger/LoggerFactory.o \
	bmpg_uncc_edu/util/xml/xmlParser.o \
	bmpg_uncc_edu/chemistry/PDBProtein.o \
	bmpg_uncc_edu/fast/FASTInitializer.o \
	AtomCoordinates.o \
	HashGridVolume.o \
	HoshenKopelman.o \
	RotationStatistics.o \
	SpringModel.o \
	main.o 

	
CalculateVolume: $(OBJS)
	$(CXX) $(CXXFLAGS) -o CalculateVolume  $(OBJS)
    
clean:
	rm -f bmpg_uncc_edu/algorithms/FastGraph.o \
	bmpg_uncc_edu/algorithms/GraphTraits.o \
	bmpg_uncc_edu/algorithms/HashGrid.o \
	bmpg_uncc_edu/chemistry/Atom.o \
	bmpg_uncc_edu/chemistry/Bond.o \
	bmpg_uncc_edu/chemistry/Complex.o \
	bmpg_uncc_edu/chemistry/CovalentBond.o \
	bmpg_uncc_edu/chemistry/CovalentBonder.o \
	bmpg_uncc_edu/chemistry/Element.o \
	bmpg_uncc_edu/chemistry/Molecule.o \
	bmpg_uncc_edu/chemistry/PDBAtom.o \
	bmpg_uncc_edu/chemistry/PDBProteinChain.o \
	bmpg_uncc_edu/chemistry/PDBProteinResidue.o \
	bmpg_uncc_edu/chemistry/detail/PDBChainManager.o \
	bmpg_uncc_edu/chemistry/helper/PDBHelper.o \
	bmpg_uncc_edu/chemistry/library/AminoAcidLibrary.o \
	bmpg_uncc_edu/chemistry/library/BondDistanceTable.o \
	bmpg_uncc_edu/chemistry/library/PeriodicTable.o \
	bmpg_uncc_edu/chemistry/library/loader/AALibLoader.o \
	bmpg_uncc_edu/fast/Defs.o \
	bmpg_uncc_edu/fast/ParameterFile.o \
	bmpg_uncc_edu/fast/input/ResourceBundle.o \
	bmpg_uncc_edu/math/MathConstants.o \
	bmpg_uncc_edu/math/R3Vector.o \
	bmpg_uncc_edu/math/Vector.o \
	bmpg_uncc_edu/util/Exception.o \
	bmpg_uncc_edu/util/FlatFileReader.o \
	bmpg_uncc_edu/util/Properties.o \
	bmpg_uncc_edu/util/SequentialIdGenerator.o \
	bmpg_uncc_edu/util/StringUtil.o \
	bmpg_uncc_edu/util/logger/DefaultLogger.o \
	bmpg_uncc_edu/util/logger/Logger.o \
	bmpg_uncc_edu/util/logger/LoggerFactory.o \
	bmpg_uncc_edu/util/xml/xmlParser.o \
	bmpg_uncc_edu/chemistry/PDBProtein.o \
	bmpg_uncc_edu/fast/FASTInitializer.o \
	AtomCoordinates.o \
	HashGridVolume.o \
	HoshenKopelman.o \
	RotationStatistics.o \
	SpringModel.o \
	main.o 

	
# smith_problem_scripts

gapHypothesis.g:
	global variables:
		complexEquivalent := NewDictionary( [], true );
		realIrreducibles := [];
		dimensionsRealModules := [];
		realIrrOfDim := [];
		numberRealIrrOfDim := [];
		realIrrNr := NewDictionary( [], true );
		realIrrNrReversed := NewDictionary( 1, true );
		realModulesGivenDimension := [];
		dimensionTable := [];
		dimensionValues := [];
	functions:
		frobeniusSchurIndicator := function( chi, G )
		realIrr := function( G )
		fixedPointDimensionIrr := function( ir, H, G )
		fixedPointDimensionRealModule := function( realModule, H, G )
		restrictedTuples := function( restrictions )
		realModulesOfDimension := function( dim )
		verifyGapHypothesisModule := function( U, W, G )
		dimensionRealModule := function( realModule )
		minimalGapModule := function( smithSetElement, G, dimensionUpperBound )
		minimalGapModuleGroup := function( G )
	
genPORed.g:
	global variables:
		complexEquivalent := NewDictionary( [], true );
		realIrreducibles := [];
		dimensionsRealModules := [];
		realIrrOfDim := [];
		numberRealIrrOfDim := [];
		realIrrNr := NewDictionary( [], true );
		realIrrNrReversed := NewDictionary( 1, true );
		realModulesGivenDimension := [];
		dimensionTable := [];
		dimensionValues := [];
	functions:
		frobeniusSchurIndicator := function( chi, G )
		realIrr := function( G )
		realModuleCharacters := function( realModule, G )
		restrictedCharacter := function( realModule, G )
		generatorsPOred := function( G )
		
smithSetPORed.g:
	global variables:
		complexEquivalent := NewDictionary( [], true );
		realIrreducibles := [];
		dimensionsRealModules := [];
		realIrrOfDim := [];
		numberRealIrrOfDim := [];
		realIrrNr := NewDictionary( [], true );
		realIrrNrReversed := NewDictionary( 1, true );
		realModulesGivenDimension := [];
		dimensionTable := [];
		dimensionValues := [];
	functions:
		frobeniusSchurIndicator := function( chi, G )
		realIrr := function( G )
		fixedPointDimensionIrr := function( ir, H, G )
		fixedPointDimensionRealModule := function( realModule, H, G )
		hasPORedAsSmithSet := function( G ) # TO VERIFY
		
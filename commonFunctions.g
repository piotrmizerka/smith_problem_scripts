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

frobeniusSchurIndicator := function( chi, G )
	local result, cl, repr;
	result := 0;
	for cl in ConjugacyClasses( G ) do
		repr := Representative( cl );
		result := result+Size( cl )*( repr*repr )^chi;
	od;
	result := result/Order( G );
	return result;
end;

realIrr := function( G )
	local irr, ir, row, cl, ind, complexIrr, complexirr, check, n, considered, i, trivialModule;
	irr := Irr( G );
	complexIrr := [];
	complexEquivalent := NewDictionary( [], true );
	realIrreducibles := [];
	dimensionsRealModules := [];
	considered := [];
	realIrrOfDim := [];
	numberRealIrrOfDim := [];
	realIrrNr := NewDictionary( [], true );
	realIrrNrReversed := NewDictionary( 1, true );
	i := 1;
	trivialModule := [];
	for cl in ConjugacyClasses( G ) do
		Add( trivialModule, 1 );
	od;
	for ir in irr do
		row := [];
		ind := frobeniusSchurIndicator( ir, G );
		if ind = 1 then
			for cl in ConjugacyClasses( G ) do
				Add( row, Representative( cl )^ir );
			od;
			if row <> trivialModule then
				Add( realIrreducibles, row );
				AddDictionary( complexEquivalent, row, ir );
				realIrrOfDim[row[1]] := [];
				considered[row[1]] := false;
				numberRealIrrOfDim[row[1]] := 0;
				AddDictionary( realIrrNr, row, i );
				AddDictionary( realIrrNrReversed, i, row );
				i := i+1;
			fi;
		elif ind = -1 then
			for cl in ConjugacyClasses( G ) do
				Add( row, 2*RealPart( Representative( cl )^ir ) );
			od;
			Add( realIrreducibles, row );
			AddDictionary( complexEquivalent, row, ir );
			realIrrOfDim[row[1]] := [];
			considered[row[1]] := false;
			numberRealIrrOfDim[row[1]] := 0;
			AddDictionary( realIrrNr, row, i );
			AddDictionary( realIrrNrReversed, i, row );
			i := i+1;
		else
			for cl in ConjugacyClasses( G ) do
				Add( row, 2*RealPart( Representative( cl )^ir ) );
			od;
			check := true;
			for complexirr in complexIrr do
				if complexirr = row then
					check := false;
					break;
				fi;
			od;
			if check = true then
				Add( realIrreducibles, row );
				Add( complexIrr, row );
				AddDictionary( complexEquivalent, row, ir );
				realIrrOfDim[row[1]] := [];
				considered[row[1]] := false;
				numberRealIrrOfDim[row[1]] := 0;
				AddDictionary( realIrrNr, row, i );
				AddDictionary( realIrrNrReversed, i, row );
				i := i+1;
			fi;
		fi;
	od;
	for ir in realIrreducibles do
		Add( realIrrOfDim[ir[1]], ir );
		numberRealIrrOfDim[ir[1]] := numberRealIrrOfDim[ir[1]]+1;
		if considered[ir[1]] = false then
			Add( dimensionsRealModules, ir[1] );
			considered[ir[1]] := true;
		fi;
	od;
end;
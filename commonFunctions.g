# File containing functions used by at least 2 different files.

# To load this file separately, paste this line (change the directory if necessary),
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "commonFunctions.g" ) );

# Global variables
complexEquivalent := NewDictionary( [], true ); # complex irreducible character corresponding to a given real irreducible
complexIrreducibles := []; # characters of complex irreducible representations
realIrreducibles := []; # characters of real irreducible representations
dimensionsRealModules := []; # dimensions of real irreducible representations
dimensionsRealModulesNoRepetitions := []; # dimensions of real irreducible representations without repetitions
realIrrOfDim := []; # real irreducible characters of a given dimenension
numberRealIrrOfDim := []; # number of real irreducible representations of a given dimension
realIrrNr := NewDictionary( [], true ); # idies of real irreducible representations
realIrrNrReversed := NewDictionary( 1, true ); # real irreducible characters for a given id

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
	complexIrreducibles := irr;
	complexIrr := [];
	complexEquivalent := NewDictionary( [], true );
	realIrreducibles := [];
	dimensionsRealModules := [];
	dimensionsRealModulesNoRepetitions := [];
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
	i := 1;
	for cl in ConjugacyClasses( G ) do
		if Order( Representative( cl ) ) = 1 then
			break;
		fi;
		i := i+1;
	od;
	for ir in realIrreducibles do
		Add( realIrrOfDim[ir[i]], ir );
		numberRealIrrOfDim[ir[i]] := numberRealIrrOfDim[ir[i]]+1;
		Add( dimensionsRealModules, ir[i] );
		if considered[ir[i]] = false then
			Add( dimensionsRealModulesNoRepetitions, ir[i] );
			considered[ir[i]] := true;
		fi;
	od;
end;

fixedPointDimensionIrr := function( ir, H, G )
	local result, h;
	result := 0;
	for h in H do
		result := result+h^ir;
	od;
	result := result/Order( H );
	if frobeniusSchurIndicator( ir, G ) <> 1 then
		result := result*2;
	fi;
	return result;
end;

fixedPointDimensionRealModule := function( realModule, H, G )
	local result, irrComponent;
	result := 0;
	for irrComponent in realModule do
		result := result+fixedPointDimensionIrr( LookupDictionary( complexEquivalent, irrComponent[1] ), H, G )*irrComponent[2];
	od;
	return result;
end;

# Requires realIrr( G ) to be called earlier
tableFixedPointDimension := function( G )
	local H, row, ir, temp, temp2, text, i, result;
	row := [];
	result := [];
	i := 1;
	for H in ConjugacyClassesSubgroups( G ) do
		Add( row, i );
		i := i+1;
	od;
	Add( result, row );
	for ir in realIrreducibles do
		row := [];
		temp := [];
		temp2 := [];
		Add( temp, ir );
		Add( temp, 1 );
		Add( temp2, temp );
		for H in ConjugacyClassesSubgroups( G ) do
			Add( row, fixedPointDimensionRealModule( temp2, Representative( H ), G ) );
		od;
		Add( result, row );
	od;
	Display( result );
	# Problems with StructureDescription - problems with SmallGroup library installation
end;

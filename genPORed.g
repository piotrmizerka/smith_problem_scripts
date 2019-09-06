# To load this file separately, paste this line (change the directory if necessary),
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "genPORed.g" ) );

# Requires "commonFunctions.g" to be called earlier. It can be called with the following command,
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "commonFunctions.g" ) );

# Computes the characters of a given real G-module.
# Requires realIrr( G ) to be called earlier
realModuleCharacters := function( realModule, G )
	local irrComponent, result, temp, cl, g, i;
	result := [];
	for i in [1..Size( ConjugacyClasses( G ) )] do
		temp := 0;
		for irrComponent in realModule do
			temp := temp+irrComponent[1][i]*irrComponent[2];
		od;
		Add( result, temp );
	od;
	return result;
end;

# Computes characters of prime power order elements of a group G for a given real G-module.
# Requires realIrr( G ) to be called earlier
restrictedCharacter := function( realModule, G ) # character restricted to elts of appropriate prime power order
	local cl, g, result, i, component, characters;
	result := [];
	i := 1;
	characters := realModuleCharacters( realModule, G );
	for cl in ConjugacyClasses( G ) do
		g := Representative( cl );
		if IsPrimePowerInt( Order( g ) ) or Order( g ) = 1 then
			Add( result, characters[i] );
		fi;
		i := i+1;
	od;
	return result;
end;

# Computes generators of the reduced primary group of G in the basis given by real irreducibles
# computed by realIrr( G ). The generators are over the real numbers, i.e. we treat the basis of
# the reduced primary group as the basis over a vector space over real numbers.
# Requires realIrr( G ) to be called earlier
generatorsPORed := function( G )
	local A, temp, irr, irrx, irrxx;
	A := [];
	realIrr( G );
	for irr in realIrreducibles do
		irrx := [];
		Add( irrx, irr );
		Add( irrx, 1 );
		irrxx := [];
		Add( irrxx, irrx );
		Add( A, restrictedCharacter( irrxx, G ) );
	od;
	return NullspaceMat( A );
	#return NullspaceMat( TransposedMat( A ) );
end;

# MAIN FUNCTION TO CALL
# Multiplies the generators of the reduced primary group of G so that they have integer coefficients
# (note that this does not imply that the so obtained generators (over R) are the generators of the
# reduced primary group - that is over the integers.
# Saves such integer "generators"-elements of the reduced primary group of G in a specified file.
# Saves dimensions of real irreducible representations in the same file.
# Requires realIrr( G ) to be called earlier.
saveGeneratorsPORed := function( G, savePath )
	local generators, integerGenerators, minFraction, generator,
				currentFraction, coefficient, i, integerGenerator, j;
	generators := generatorsPORed( G );
	integerGenerators := [];
	for generator in generators do
		minFraction := 1;
		for coefficient in generator do
			currentFraction := coefficient-Int( coefficient );
			if minFraction > currentFraction and currentFraction > 0 then
				minFraction := currentFraction;
			fi;
		od;
		integerGenerator := [];
		for i in [1..Size( generator )] do
				Add( integerGenerator, Int( generator[i]/minFraction ) );
		od;
		Add( integerGenerators, integerGenerator );
	od;
	PrintTo( savePath, "genaratorVerticesFromFile = " );
	AppendTo( savePath, "[ " );
	for i in [1..Size( integerGenerators )] do
		AppendTo( savePath, "( ");
		for j in [1..Size( integerGenerators[i] )-1] do
			AppendTo( savePath, integerGenerators[i][j] );
			AppendTo( savePath, ", " );
		od;
		AppendTo( savePath, integerGenerators[i][Size(integerGenerators[i])] );
		AppendTo( savePath, " )" );
		if i < Size( integerGenerators ) then
			AppendTo( savePath, ", " );
		fi;
	od;
	AppendTo( savePath, " ]\n" );
	AppendTo( savePath, "realIrreduciblesDimensionsFromFile = ( " );
	for i in [1..Size( dimensionsRealModules )] do
		AppendTo( savePath, dimensionsRealModules[i] );
		if i < Size( dimensionsRealModules ) then
			AppendTo( savePath, ", " );
		fi;
	od;
	AppendTo( savePath, " )" );
end;

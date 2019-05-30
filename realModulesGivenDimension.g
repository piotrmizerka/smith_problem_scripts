# Functions responsible for determining all the real modules of a group G of a given dimension.

# To load this file paste this line (change the directory if necessary),
# Read( Filename( DirectoryDesktop(), "realModulesGivenDimension.g" ) );

Read( Filename( DirectoryDesktop(), "commonFunctions.g" ) );

# Global variables
realModulesGivenDimension := []; # RG-modules of a given dimension

# Auxiliary function computing for a given tuple of positive integers, [a1,...,an]
# all the tuples [b1,...,bn] of positive integers with bi<=ai for i=1,...,n
restrictedTuples := function( restrictions )
	local result, tuple, finalTuple, res, coord, value, i, temp, tuple2;
	finalTuple := [];
	tuple := [];
	coord := Size( restrictions );
	value := 1;
	result := [];
	for res in restrictions do
		Add( finalTuple, res );
		Add( tuple, 1 );
	od;
	while tuple <> finalTuple do
		tuple2 := [];
		for i in [1..Size(restrictions)] do
			Add( tuple2, tuple[i] );
		od;
		Add( result, tuple2 );
		while tuple[coord]+1 > restrictions[coord] do
			coord := coord-1;
		od;
		tuple[coord] := tuple[coord]+1;
		for i in [(coord+1)..Size(restrictions)] do
			tuple[i] := 1;
		od;
		coord := Size( restrictions );
	od;
	Add( result, finalTuple );
	return result;
end;

# Main function.
# Computes all the real modules of a given dimension.
# It saves the result in the global variable "realModulesGivenDimension".
# Requires realIrr( G ) to be called earlier.
realModulesOfDimension := function( dim )
	local restrictedPartitions, numberInPartitionDim, n, partition, summand, summands, 
		  restrictions, set, ir, resTup, restup, unTup, uT, result, realModule, 
		  multiplicities, tempModule, temp, temp2, i;
	numberInPartitionDim := [];
	restrictedPartitions := RestrictedPartitions( dim, dimensionsRealModules );
	realModulesGivenDimension := [];
	result := [];
	for partition in restrictedPartitions do
		for n in dimensionsRealModules do
			numberInPartitionDim[n] := 0;
		od;
		summands := [];
		for summand in partition do
			numberInPartitionDim[summand] := numberInPartitionDim[summand]+1;
				if ( summand in summands ) = false then
					Add( summands, summand );
				fi;
		od;
		restrictions := [];
		unTup := [];
		for summand in summands do
			set := [];
			for ir in realIrrOfDim[summand] do
				Add( set, ir );
			od;
			uT := [];
			uT := UnorderedTuples( set, numberInPartitionDim[summand] );
			unTup[summand] := uT;
			Add( restrictions, Size( unTup[summand] ) );
		od;
		resTup := restrictedTuples( restrictions );
		for restup in resTup do
			multiplicities := [];
			realModule := [];
			for ir in realIrreducibles do
				Add( multiplicities, 0 );
			od;
			for i in [1..Size(summands)] do
				for tempModule in unTup[summands[i]][restup[i]] do
					temp := LookupDictionary( realIrrNr, tempModule );
					multiplicities[temp] := multiplicities[temp]+1;
				od;
			od;
			for i in [1..Size(realIrreducibles)] do
				if multiplicities[i] > 0 then
					temp2 := [];
					Add( temp2, LookupDictionary( realIrrNrReversed, i ) );
					Add( temp2, multiplicities[i] );
					Add( realModule, temp2 );
				fi;
			od;
			Add( result, realModule );
		od;
	od;
	realModulesGivenDimension := result;
end;
# Functions responsible for finding an RG-modules V satisfying the following conditions,
# (1) dimV^P >= 2dimV^H for any p-subgroup P of G and P<H<=G (weak gap hypothesis),
# (2) dimV^P >= 5 and dimV^H >= 2 for any p-subgroup P of G and pseudocyclic subgroup H of G,
# (3) every pseudocyclic subgroup occurs as an isotropy subgroup of the action of G 
#	  on the considered RG-module V,
# (4) V^L = {0} for every large subgroup of G,
# (5) V is P-orientable.
# For the definitions of the above concepts, see [1].

# To load this file paste this line (change the directory if necessary),
# Read( Filename( DirectoryDesktop(), "sufficientConditions.g" ) );

Read( Filename( DirectoryDesktop(), "realModulesGivenDimension.g" ) );

# Global variables
condition1SubgroupPairs := []; # pairs of subgroups of a group G of the form (P,H) where P is a p-group and P<H
pSubgroups := []; # p-subgroups of G
pseudocyclicSubgroups :=[]; # pseudocyclic subgroups of G
condition3SubgroupPairs := []; # pairs of subgroups of a group G of the form (H,K) where H is pseudocyclic and H<K
smallestNormalSubgroupsQuotientPGroup := NewDictionary( 1, true ); # OpG subgroups for a group G - for the definition, see [1]
largeSubgroups := []; # large subgroups of G
elementsOfPrimePowerOrder := []; # elements of G of prime power order
complexIrreducibleRepresentations := []; # complex irreducible representations affording complex irreducible characters
# Remark: all the subgroup tuples are determined up to conjugacy since it does not affect fixed point dimensions.

# Condition (1) function.
# Computes subgroup pairs (P,H) where P<H<=G and P is a p-group.
# Saves the result in the "condition1SubgroupPairs".
determineCondition1SubgroupPairs := function( G )
	local clG, clH, H, P;
	condition1SubgroupPairs := [];
	for clG in ConjugacyClassesSubgroups( G ) do
		H := Representative( clG );
		for clH in ConjugacyClassesSubgroups( H ) do
			P := Representative( clH );
			if IsPGroup( P ) and Order( P ) < Order( H ) then
				Add( condition1SubgroupPairs, [P,H] );
			fi;
		od;
	od;
end;

# Condition (1) main function. # TO VERIFY
# Checks the weak gap hypothesis (condition (1)) for an RG-module U+W.
# Requires realIrr( G ) and determineCondition1SubgroupPairs( G ) to be called earlier.
verifyGapHypothesisModule := function( U, W, G )
	local PHPair, P, H;
	for PHPair in condition1SubgroupPairs do
		P := PHPair[1];
		H := PHPair[2];
		if fixedPointDimensionRealModule( U, P, G )+fixedPointDimensionRealModule( W, P, G ) < 
		   2*(fixedPointDimensionRealModule( U, H, G )+fixedPointDimensionRealModule( W, H, G )) then
			return false;
		fi;
	od;
	return true;
end;

# Condition (2) auxilary function.
# Computes the quotient G/N for a given normal subgroup N of G.
quotientGroup := function( G, N )
	return Image( NaturalHomomorphismByNormalSubgroup( G, N ) );
end;

# Conditon (2) and (3) auxiliary function.
# Checks if a group G is pseudocyclic.
isPseudocyclic := function( G )
	local N;
	for N in NormalSubgroups( G ) do
		if (IsPGroup( N ) = true) and (IsCyclic( quotientGroup( G, N ) ) = true) then
			return true;
		fi;
	od;
	return false;
end;

# Condition (2) and (3) function.
# Determines all pseudocylcic subgroups of G up to conjugacy.
# Saves the result in the global variable "pseudocyclicSubgroups".
determinePseudocyclicSubgroups := function( G )
	local cl, H;
	pseudocyclicSubgroups := [];
	for cl in ConjugacyClassesSubgroups( G ) do
		H := Representative( cl );
		if isPseudocyclic( H ) then
			Add( pseudocyclicSubgroups, H );
		fi;
	od;
end;

# Condition (2) function.
# Determines all the p-subgroups of G up to conjugacy.
# Saves the result in the global variable "pSubgroups".
determinePSubgroups := function( G )
	local cl, H;
	pSubgroups := [];
	for cl in ConjugacyClassesSubgroups( G ) do
		H := Representative( cl );
		if IsPGroup( H ) then
			Add( pSubgroups, H );
		fi;
	od;
end;

# Condition (2) main function.
# Checks the condition (2) for an RG-module U+W # TO VERIFY
# Requires realIrr( G ), determinePseudocyclicSubgroups( G ) and determinePSubgroups( G ) to be called earlier.
verifyCondition2 := function( U, W, G )
	local P, H;
	for P in pSubgroups do
		if fixedPointDimensionRealModule( U, P, G )+fixedPointDimensionRealModule( W, P, G ) < 5 then
			#Display( U );
			#Display( W );
			#Print( fixedPointDimensionRealModule( U, P, G )+fixedPointDimensionRealModule( W, P, G ), "a\n" );
			return false;
		fi;
	od;
	for H in pseudocyclicSubgroups do
		if fixedPointDimensionRealModule( U, H, G )+fixedPointDimensionRealModule( W, H, G ) < 2 then
			#Display( U );
			#Display( W );
			#Print( fixedPointDimensionRealModule( U, H, G )+fixedPointDimensionRealModule( W, H, G ), "b\n" );
			return false;
		fi;
	od;
	return true;
end;

# Condition (3) function.
# Determines pairs of subgroups of G (up to conjugacy) of the form [H,K] where H i pseudocyclic and H<K.
# Saves the result in the global variable "condition3SubgroupPairs".
determineCondition3SubgroupPairs := function( G )
	local clG, clK, H, K;
	condition3SubgroupPairs := [];
	for clG in ConjugacyClassesSubgroups( G ) do
		K := Representative( clG );
		for clK in ConjugacyClassesSubgroups( K) do
			H := Representative( clK );
			if isPseudocyclic( H ) and Order( H ) < Order( K ) then
				Add( condition3SubgroupPairs, [H,K] );
			fi;
		od;
	od;
end;

# Condition (3) main function. # TO VERIFY
# Checks the condition (3) for an RG-module U+W. 
# Requires realIrr( G ) and determineCondition3SubgroupPairs( G ) to be called earlier.
verifyCondition3 := function( U, W, G )
	local HKPair, H, K;
	for HKPair in condition3SubgroupPairs do
		H := HKPair[1];
		K := HKPair[2];
		if fixedPointDimensionRealModule( U, H, G )+fixedPointDimensionRealModule( W, H, G ) = 
		   fixedPointDimensionRealModule( U, K, G )+fixedPointDimensionRealModule( W, K, G ) then
			return false;
		fi;
	od;
	return true;
end;

# Condition (4) function.
# Computes smallest normal subgroups of such that the quotients are p-groups for primes p dividing |G|.
# Saves the result in the global variable "smallestNormalSubgroupsQuotientPGroup".
determineSmallestNormalSubgroupsQuotientPGroup := function( G )
	local N, quotient, smallestOrders, quotientOrder, p, primeDivisors, subgroupOrder;
	smallestNormalSubgroupsQuotientPGroup := NewDictionary( 1, true ); 
	smallestOrders := NewDictionary( 1, true );
	for p in PrimeDivisors( Order( G ) ) do
		AddDictionary( smallestOrders, p, 1000000 );
	od;
	for N in NormalSubgroups( G ) do
		quotient := quotientGroup( G, N );
		quotientOrder := Order( quotient );
		primeDivisors := PrimeDivisors( quotientOrder );
		subgroupOrder := Order( N );
		if Size( primeDivisors ) > 0 then
			p := primeDivisors[1];
		else
			p := 1;
		fi;
		if IsPGroup( quotient ) and LookupDictionary( smallestOrders, p ) > subgroupOrder  then
			AddDictionary( smallestNormalSubgroupsQuotientPGroup, p, N );
			AddDictionary( smallestOrders, p, subgroupOrder );
		fi;
	od;
end;

# Condition (4) function.
# Checks if a subgroup L of a group G is its large subgroup
# (for the definition of a large subgroup, see [1]).
# Requires "determineSmallestNormalSubgroupsQuotientPGroup( G )" to be called earlier.
isLargeSubgroup := function( L, G )
	local p, OpG;
	for p in PrimeDivisors( Order( G ) ) do
		OpG := LookupDictionary( smallestNormalSubgroupsQuotientPGroup, p );
		if IsSubgroup( L, OpG ) then
			return true;
		fi;
	od;
	return false;
end;

# Condition (4) function.
# Determines the large subgroups of a group G which are the representatives of 
# its conjgacy classes of subgroups (we restrict to conjugacy classes since the fixed point dimension is constant on them).
# Saves the result in the global variable "largeSubgroups" which is a list of subgroups of G.
determineLargeSubgroups := function( G )
	local cl, L;
	determineSmallestNormalSubgroupsQuotientPGroup( G );
	largeSubgroups := [];
	for cl in ConjugacyClassesSubgroups( G ) do
		L := Representative( cl );
		if isLargeSubgroup( L, G ) then
			Add( largeSubgroups, L );
		fi;
	od;
end;

# Condition (4) main function.
# Checks the condition (4) for for an RG-module U+W.
# Requires realIrr( G ) and determineLargeSubgroups( G ) to be called earlier.
verifyCondition4 := function( U, W, G )
	local L;
	for L in largeSubgroups do
		if fixedPointDimensionRealModule( U, L, G )+fixedPointDimensionRealModule( W, L, G ) > 0 then
			return false;
		fi;
	od;
	return true;
end;

# Condition (5) function.
# Determines elements of G of prime power order up to conjugacy.
determineElementsOfPrimePowerOrder := function( G )
	local cl, g;
	elementsOfPrimePowerOrder := [];
	for cl in ConjugacyClasses( G ) do
		g := Representative( cl );
		if IsPrimePowerInt( Order( g ) ) or Order( g ) = 1 then
			Add( elementsOfPrimePowerOrder, g );
		fi;
	od;
end;

# Condition (5) function.
# Determines complex irreducible representations affording complex irreducible characters.
determineComplexIrreducibleRepresentations := function( G )
	local irr;
	complexIrreducibleRepresentations := [];
	for irr in complexIrreducibles do
		Add( complexIrreducibleRepresentations, IrreducibleRepresentationsDixon( G, irr ) );
	od;
end;

# Condition (5) main function. # TO VERIFY
# Checks the P-orientability for an RG module U+W.
# Requires realIrr( G ), determineElementsOfPrimePowerOrder( G ) and 
# determineComplexIrreducibleRepresentations( G ) to be called earlier.
verifyPOrientabilityModule := function( U, W, G )
	local g, determinant, component, characterUW, i, coefficients, complexIrrRep;
	characterUW := [];
	for i in [1..Size( realIrreducibles )] do
		Add( characterUW, 0 );
	od;
	for component in U do
		for i in [1..Size( component[1] )] do
			characterUW[i] := characterUW[i]+component[1][i]*component[2];
		od;
	od;
	for component in W do
		for i in [1..Size( component[1] )] do
			characterUW[i] := characterUW[i]+component[1][i]*component[2];
		od;
	od;
	coefficients := SolutionMat( complexIrreducibles, characterUW );
	i := 1;
	for g in elementsOfPrimePowerOrder do
		determinant := 1;
		for complexIrrRep in complexIrreducibleRepresentations do
			determinant := determinant*DeterminantMat( Image( complexIrrRep, g ) )^coefficients[i];
		od;
		if determinant < 0 then
			return false;
		fi;
		i := i+1;
	od;
	return true;
end;

# Computes the dimension of a real module
dimensionRealModule := function( realModule )
	local result, component;
	result := 0;
	for component in realModule do
		result := result+component[1][1]*component[2];
	od;
	return result;
end;

# Initializes subgroup tuples for G
initializeSubgroupTuples := function( G )
	determineCondition1SubgroupPairs( G );
	determineCondition3SubgroupPairs( G );
	determinePSubgroups( G );
	determinePseudocyclicSubgroups( G );
	determineLargeSubgroups( G );
	determineElementsOfPrimePowerOrder( G );
	determineComplexIrreducibleRepresentations( G );
end;

# For a given smith set element U-V, computes W s.t. U+W and V+W # TO VERIFY
# satisfy the conditions (1)-(4) and dim(U+W)=dim(V+W) is minimized.
# The "smithSetElement" U-V is given in the base of RO(G) which are the real irreducible representations of G.
# U-V is read from the file "smithSetsElements.g" which should be in the Desktop directory (chagne if necessary).
# Returns the result in the list format [dim(U+W),U,V,W].
# Requires realIrr( G )and initializeSubgroupTuples( G ) to be called earlier.
minimalGapModule := function( smithSetElement, G, dimensionUpperBound )
	local U, V, W, n, coord, i;
	U := [];
	V := [];
	i := 1;
	for coord in smithSetElement do
		if coord > 0 then
			Add( U, [realIrreducibles[i],coord] );
		elif coord < 0 then
			Add( V, [realIrreducibles[i],-coord] );
		fi;
		i := i+1;
	od;
	if verifyGapHypothesisModule( U, [], G ) = true and verifyGapHypothesisModule( V, [], G ) = true and
	   verifyCondition2( U, [], G ) = true and verifyCondition2( V, [], G ) = true and
	   verifyCondition3( U, [], G ) = true and verifyCondition3( V, [], G ) = true and
	   verifyCondition4( U, [], G ) = true and verifyCondition4( V, [], G ) = true and
	   (verifyPOrientabilityModule( U, [], G ) = true or verifyPOrientabilityModule( V, [], G ) = true) then
		return [dimensionRealModule( U ),U,V,[]];
	fi;
	for n in [1..(dimensionUpperBound-dimensionRealModule( U ))] do
		#Display( n );
		realModulesOfDimension( n );
		for W in realModulesGivenDimension do
			if verifyGapHypothesisModule( U, W, G ) = true and verifyGapHypothesisModule( V, W, G ) = true and
			   verifyCondition2( U, W, G ) = true and verifyCondition2( V, W, G ) = true and
			   verifyCondition3( U, W, G ) = true and verifyCondition3( V, W, G ) = true and
			   verifyCondition4( U, W, G ) = true and verifyCondition4( V, W, G ) = true and
			   (verifyPOrientabilityModule( U, W, G ) = true or verifyPOrientabilityModule( V, W, G ) = true) then
				return [dimensionRealModule( U )+n,U,V,W];
			fi;
		od;
	od;
	return [];
end;

# Main function to call. 
# Computes RG-modules U+W and V+W which satisfy the conditions (1)-(4).
# Requires smith set elements to be already computed
# and saved in the file "smithSetElements.g" in DirectoryDesktop() (change if necessary).
minimalGapModuleGroup := function( G )
	local i, gapRealization, minGapRealization;
	realIrr( G );
	initializeSubgroupTuples( G );
	minGapRealization := [100,[],[],[]]; 
	Read( Filename( DirectoryDesktop(), "smithSetElements.g" ) );
	Print( "Number of PO red elts to consider: ", Size( smithSetElements ), "\n" );
	for i in [1..Size( smithSetElements )] do
		Print( "Id of element: ", i, "\n" );
		if smithSetElements[i][2] >= minGapRealization[1] then
			break;
		fi;
		gapRealization := minimalGapModule( smithSetElements[i][1], G, minGapRealization[1] );
		if Size( gapRealization ) > 0 then
			if gapRealization[1] < minGapRealization[1] then
				minGapRealization := gapRealization;
			fi;
		fi;
	od;
	return minGapRealization;
end;

# REFERENCES
# [1] M. Morimoto, K. Pawałowski, Smooth actions of finite Oliver groups on spheres,
#     Topology 42 (2003) 395-421.
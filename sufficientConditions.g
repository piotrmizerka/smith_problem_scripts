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
constraintsCondition1 := [];
constraintsCondition2 := [];
constraintsCondition3 := [];
constraintsCondition4 := [];
constraints := [];
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
# Checks if a proper subgroup L of a group G is its large subgroup
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
# Determines the large proper subgroups of a group G which are the representatives of 
# its conjgacy classes of subgroups (we restrict to conjugacy classes since the fixed point dimension is constant on them).
# Saves the result in the global variable "largeSubgroups" which is a list of subgroups of G.
# Requires "determineSmallestNormalSubgroupsQuotientPGroup( G )" to be called earlier.
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
	for i in [1..Size( realIrreducibles[1] )] do
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

# Determines constraints for weak gap hyopthesis for integer linear programming.
# Requires realIrr( G ) and determineCondition1SubgroupPairs( G ) to be called earlier
constraintsGapHypothesis := function( G )
	local pair, H, P, irr, constraint, dimP, dimH;
	constraintsCondition1 := [];
	for pair in condition1SubgroupPairs do
		P := pair[1];
		H := pair[2];
		constraint := [];
		for irr in realIrreducibles do
			dimP := fixedPointDimensionRealModule( [[irr,1]], P, G );
			dimH := fixedPointDimensionRealModule( [[irr,1]], H, G );
			Add( constraint, dimP-2*dimH );
		od;
		Add( constraint, ">=" );
		Add( constraint, 0 );
		Add( constraintsCondition1, constraint );
	od;
end;

# Determines constraints for the condition 2 for integer linear programming.
# Requires realIrr( G ), determinePseudocyclicSubgroups( G ) and
# determinePSubgroups( G ) to be called earlier
constraintsCondition2 := function( G )
	local H, P, irr, constraint, dim;
	constraintsCondition2 := [];
	for P in pSubgroups do
		constraint := [];
		for irr in realIrreducibles do
			dim := fixedPointDimensionRealModule( [[irr,1]], P, G );
			Add( constraint, dim );
		od;
		Add( constraint, ">=" );
		Add( constraint, 5 );
		Add( constraintsCondition2, constraint );
	od;
	for H in pseudocyclicSubgroups do
		constraint := [];
		for irr in realIrreducibles do
			dim := fixedPointDimensionRealModule( [[irr,1]], H, G );
			Add( constraint, dim );
		od;
		Add( constraint, ">=" );
		Add( constraint, 2 );
		Add( constraintsCondition2, constraint );
	od;
end;

# Determines constraints for the condition 3 for integer linear programming.
# Requires realIrr( G ) and determineCondition1SubgroupPairs( G ) to be called earlier
constraintsCondition3 := function( G )
	local pair, H, K, irr, constraint, dimH, dimK;
	constraintsCondition3 := [];
	for pair in condition3SubgroupPairs do
		H := pair[1];
		K := pair[2];
		constraint := [];
		for irr in realIrreducibles do
			dimH := fixedPointDimensionRealModule( [[irr,1]], H, G );
			dimK := fixedPointDimensionRealModule( [[irr,1]], K, G );
			Add( constraint, dimH-dimK );
		od;
		Add( constraint, ">=" );
		Add( constraint, 1 );
		Add( constraintsCondition3, constraint );
	od;
end;

# Determines constraints for the condition 4 for integer linear programming.
# Requires realIrr( G ), "determineSmallestNormalSubgroupsQuotientPGroup( G )"
# and determineLargeSubgroups( G ) to be called earlier (in that order).
constraintsCondition4 := function( G )
	local L, irr, constraint, dim;
	constraintsCondition4 := [];
	for L in largeSubgroups do
		constraint := [];
		for irr in realIrreducibles do
			dim := fixedPointDimensionRealModule( [[irr,1]], L, G );
			Add( constraint, dim );
		od;
		Add( constraint, "==" );
		Add( constraint, 0 );
		Add( constraintsCondition4, constraint );
	od;
end;

# Requires realIrr( G ) to be called earlier
determineConstraints := function( G )
	determineCondition1SubgroupPairs( G );
	determineCondition3SubgroupPairs( G );
	determinePSubgroups( G );
	determinePseudocyclicSubgroups( G );
	determineSmallestNormalSubgroupsQuotientPGroup( G );
	determineLargeSubgroups( G );
	determineElementsOfPrimePowerOrder( G );
	determineComplexIrreducibleRepresentations( G );
	constraints := [];
	constraintsGapHypothesis( G );
	constraintsCondition2( G );
	constraintsCondition3( G );
	constraintsCondition4( G );
	Add( constraints, constraintsCondition1 );
	Add( constraints, constraintsCondition2 );
	Add( constraints, constraintsCondition3 );
	Add( constraints, constraintsCondition4 );
end;

# TO VERIFY
# Requires realIrr( G ) and determineConstraints( G ) to be called earlier.
saveConstraintsAsSAGEFile := function( savePath, smithSetElement )
	local constraint, i, j, U, V, Ux, Vx, P, H, dimP, dimH, rhs;
	U := [];
	V := [];
	Ux := [];
	Vx := [];
	for i in [1..Size( smithSetElement )] do
		if smithSetElement[i] > 0 then
			U[i] := smithSetElement[i];
			V[i] := 0;
			Add( Ux, [realIrreducibles[i],smithSetElement[i]] );
		else
			U[i] := 0;
			V[i] := -smithSetElement[i];
			Add( Vx, [realIrreducibles[i],-smithSetElement[i]] );
		fi;
	od;
	PrintTo( savePath, "equationSystem = MixedIntegerLinearProgram( maximization=False, solver = \"GLPK\" )\n" );
	AppendTo( savePath, "x = equationSystem.new_variable( integer = True, nonnegative = True )\n\n\n" );
	AppendTo( savePath, "#Constraints from gap hypothesis\n" );
	for j in [1..Size( constraintsCondition1 )] do
		appendConstraint1( savePath, j, Ux );
		appendConstraint1( savePath, j, Vx );
	od;
	AppendTo( savePath, "\n#Constraints from condition 2\n" );
	for j in [1..Size( pSubgroups )] do
		appendConstraint2( savePath, j, (true), j, Ux );
		appendConstraint2( savePath, j, (true), j, Vx );		
	od;
	for j in [1..Size( pseudocyclicSubgroups )] do
		appendConstraint2( savePath, j+Size( pSubgroups ), false, j, Ux );
		appendConstraint2( savePath, j+Size( pSubgroups ), false, j, Vx );		
	od;
	AppendTo( savePath, "\n#Constraints from condition 3\n" );
	for j in [1..Size( constraintsCondition3 )] do
		appendConstraint3( savePath, j, Ux );
		appendConstraint3( savePath, j, Vx );
	od;
	AppendTo( savePath, "\n#Constraints from condition 4\n" );
	for j in [1..Size( constraintsCondition4 )] do
		appendConstraint4( savePath, j, Ux );
		appendConstraint4( savePath, j, Vx );
	od;
	AppendTo( savePath, "\n\nequationSystem.set_objective( " );
	for i in [1..Size( realIrreducibles )-1] do
		AppendTo( savePath, realIrreducibles[i][1] );
		AppendTo( savePath, "*x[" );
		AppendTo( savePath, i-1 );
		AppendTo( savePath, "]+" );
	od;
	AppendTo( savePath, realIrreducibles[Size( realIrreducibles )][1] );
	AppendTo( savePath, "*x[" );
	AppendTo( savePath, Size( realIrreducibles )-1 );
	AppendTo( savePath, "] )\n\nequationSystemshow()\n" );
	AppendTo( savePath, "\n\nprint( \'Objective Value: {}\'.format( equationSystem.solve() ) )\n" );
	AppendTo( savePath, "for i, v in sorted( equationSystem.get_values( x ).items()):\n" );
	AppendTo( savePath, "\tprint( \'w_\%s = \%s\' \% (i, int( round( v ) )) )" );
end;

appendConstraint1 := function( savePath, constraintId, realModule )
	local i, P, H, dimP, dimH, rhs;
	AppendTo( savePath, "equationSystem.add_constraint( " );
	for i in [1..Size( constraintsCondition1[constraintId] )-3] do
		AppendTo( savePath, "(" );
		AppendTo( savePath, constraintsCondition1[constraintId][i] );
		AppendTo( savePath, ")" );
		AppendTo( savePath, "*x[" );
		AppendTo( savePath, i-1 );
		AppendTo( savePath, "]+" );
	od;
	AppendTo( savePath, "(" );
	AppendTo( savePath, constraintsCondition1[constraintId][Size( constraintsCondition1[constraintId] )-2] );
	AppendTo( savePath, ")" );
	AppendTo( savePath, "*x[" );
	AppendTo( savePath, Size( constraintsCondition1[constraintId] )-3 );
	AppendTo( savePath, "] " );
	AppendTo( savePath, constraintsCondition1[constraintId][Size( constraintsCondition1[constraintId] )-1] );
	AppendTo( savePath, " " );
	P := condition1SubgroupPairs[constraintId][1];
	H := condition1SubgroupPairs[constraintId][2];
	dimP := fixedPointDimensionRealModule( realModule, P, G );
	dimH := fixedPointDimensionRealModule( realModule, H, G );
	rhs := -2*dimP+dimH+constraintsCondition1[constraintId][Size( constraintsCondition1[constraintId] )];
	AppendTo( savePath, rhs );
	AppendTo( savePath, " )\n" );
end;

appendConstraint2 := function( savePath, constraintId, check, subgroupId, realModule )
	local i, P, H, dimP, dimH, rhs;
	AppendTo( savePath, "equationSystem.add_constraint( " );
	for i in [1..Size( constraintsCondition2[constraintId] )-3] do
		AppendTo( savePath, "(" );
		AppendTo( savePath, constraintsCondition2[constraintId][i] );
		AppendTo( savePath, ")" );
		AppendTo( savePath, "*x[" );
		AppendTo( savePath, i-1 );
		AppendTo( savePath, "]+" );
	od;
	AppendTo( savePath, "(" );
	AppendTo( savePath, constraintsCondition2[constraintId][Size( constraintsCondition2[constraintId] )-2] );
	AppendTo( savePath, ")" );
	AppendTo( savePath, "*x[" );
	AppendTo( savePath, Size( constraintsCondition2[constraintId] )-3 );
	AppendTo( savePath, "] " );
	AppendTo( savePath, constraintsCondition2[constraintId][Size( constraintsCondition2[constraintId] )-1] );
	AppendTo( savePath, " " );
	if check = true then
		P := pSubgroups[subgroupId];
		dimP := fixedPointDimensionRealModule( realModule, P, G );
		rhs := 5-dimP;
	else
		H := pseudocyclicSubgroups[subgroupId];
		dimH := fixedPointDimensionRealModule( realModule, H, G );
		rhs := 2-dimH;
	fi;
	AppendTo( savePath, rhs );
	AppendTo( savePath, " )\n" );
end;

appendConstraint3 := function( savePath, constraintId, realModule )
	local i, H, K, dimH, dimK, rhs;
	AppendTo( savePath, "equationSystem.add_constraint( " );
	for i in [1..Size( constraintsCondition3[constraintId] )-3] do
		AppendTo( savePath, "(" );
		AppendTo( savePath, constraintsCondition3[constraintId][i] );
		AppendTo( savePath, ")" );
		AppendTo( savePath, "*x[" );
		AppendTo( savePath, i-1 );
		AppendTo( savePath, "]+" );
	od;
	AppendTo( savePath, "(" );
	AppendTo( savePath, constraintsCondition3[constraintId][Size( constraintsCondition3[constraintId] )-2] );
	AppendTo( savePath, ")" );
	AppendTo( savePath, "*x[" );
	AppendTo( savePath, Size( constraintsCondition3[constraintId] )-3 );
	AppendTo( savePath, "] " );
	AppendTo( savePath, constraintsCondition3[constraintId][Size( constraintsCondition3[constraintId] )-1] );
	AppendTo( savePath, " " );
	H := condition3SubgroupPairs[constraintId][1];
	K := condition3SubgroupPairs[constraintId][2];
	dimH := fixedPointDimensionRealModule( realModule, H, G );
	dimK := fixedPointDimensionRealModule( realModule, K, G );
	rhs := 1-dimH+dimK;
	AppendTo( savePath, rhs );
	AppendTo( savePath, " )\n" );
end;

appendConstraint4 := function( savePath, constraintId, realModule )
	local i, L, dimL, rhs;
	AppendTo( savePath, "equationSystem.add_constraint( " );
	for i in [1..Size( constraintsCondition4[constraintId] )-3] do
		AppendTo( savePath, "(" );
		AppendTo( savePath, constraintsCondition4[constraintId][i] );
		AppendTo( savePath, ")" );
		AppendTo( savePath, "*x[" );
		AppendTo( savePath, i-1 );
		AppendTo( savePath, "]+" );
	od;
	AppendTo( savePath, "(" );
	AppendTo( savePath, constraintsCondition4[constraintId][Size( constraintsCondition4[constraintId] )-2] );
	AppendTo( savePath, ")" );
	AppendTo( savePath, "*x[" );
	AppendTo( savePath, Size( constraintsCondition4[constraintId] )-3 );
	AppendTo( savePath, "] " );
	AppendTo( savePath, constraintsCondition4[constraintId][Size( constraintsCondition4[constraintId] )-1] );
	AppendTo( savePath, " " );
	L := largeSubgroups[constraintId];
	dimL := fixedPointDimensionRealModule( realModule, L, G );
	rhs := -dimL;
	AppendTo( savePath, rhs );
	AppendTo( savePath, " )\n" );
end;

# EXAMPLE
G := SL( 2, 5 );
realIrr( G );
determineConstraints( G );
saveConstraintsAsSAGEFile( "C:\\Users\\Piotrek\\Desktop\\test.sage", [-1,1,2,-2,0,0,0,0] );

# REFERENCES
# [1] M. Morimoto, K. Pawałowski, Smooth actions of finite Oliver groups on spheres,
#     Topology 42 (2003) 395-421.
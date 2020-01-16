# To load this file separately, paste this line (change the directory if necessary),
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "induction.g" ) );

Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "commonFunctions.g" ) );

isComplexInductionHomomorphismInjective := function( H, G )
	local h, sizeClH, sizeClGIntersectionH;
	for h in H do
		sizeClH := Size( ConjugacyClass( H, h ) );
		sizeClGIntersectionH := Size( Intersection( ConjugacyClass( G, h ), H ) );
		if sizeClH < sizeClGIntersectionH then
			return false;
		fi;
	od;
	return true;
end;

normalSubgroupsWithComplexInductionMonomorphism := function( G )
	local N;
	for N in NormalSubgroups( G ) do
		if isComplexInductionHomomorphismInjective( N, G ) = true then
			Print( IdGroup( N ), " = ", StructureDescription( N ), "\n" );
		fi;
	od;
end;

subgroupsWithComplexInductionMonomorphism := function( G )
	local H;
	for H in AllSubgroups( G ) do
		if isComplexInductionHomomorphismInjective( H, G ) = true then
			Print( IdGroup( H ), " = ", StructureDescription( H ), "\n" );
		fi;
	od;
end;

isRealInductionHomomorphismInjective := function( H, G )
	local h, sizeRealClH, sizeRealClGIntersectionH, int1, int2;
	for h in H do
		sizeRealClH := Size( ConjugacyClass( H, h ) )+Size( ConjugacyClass( H, h^(-1) ) )-
					   Size( Intersection( ConjugacyClass( H, h ), ConjugacyClass( H, h^(-1) ) ) );
		int1  := Intersection( ConjugacyClass( G, h ), H );
		int2  := Intersection( ConjugacyClass( G, h^(-1) ), H );
		sizeRealClGIntersectionH := Size( int1 )+Size( int2 )-Size( Intersection( int1, int2 ) );
		if sizeRealClH < sizeRealClGIntersectionH then
			return false;
		fi;
	od;
	return true;
end;

normalSubgroupsWithRealInductionMonomorphism := function( G )
	local N;
	for N in NormalSubgroups( G ) do
		if isRealInductionHomomorphismInjective( N, G ) = true then
			Print( IdGroup( N ), " = ", StructureDescription( N ), "\n" );
		fi;
	od;
end;

subgroupsWithRealInductionMonomorphism := function( G )
	local H;
	for H in AllSubgroups( G ) do
		if isRealInductionHomomorphismInjective( H, G ) = true then
			Print( IdGroup( H ), " = ", StructureDescription( H ), "\n" );
		fi;
	od;
end;

inducedCharacter := function( g, irH, H, G )
	local cosetRepresentatives, gi, result, temp;
	cosetRepresentatives := RightTransversal( G, H );
	result := 0;
	for gi in cosetRepresentatives do
		temp := gi*g*gi^(-1);
		if temp in H then
			result := result+temp^irH;
		fi;
	od;
	return result;
end;

fixedPointDimensionIrrInduced := function( irH, H, G, K )
	local result, k;
	result := 0;
	for k in K do
		result := result+inducedCharacter( k, irH, H, G );
	od;
	result := result/Order( K );
	if frobeniusSchurIndicator( irH, H ) <> 1 then
		result := result*2;
	fi;
	return result;
end;

fixedPointDimensionRealModuleInduced := function( realModuleH, H, G, K )
	local result, irrComponent;
	result := 0;
	for irrComponent in realModuleH do
		result := result+fixedPointDimensionIrrInduced( LookupDictionary( complexEquivalent, irrComponent[1] ), H, G, K )*irrComponent[2];
	od;
	return result;
end;

inducedDimensionsSubgroups := function( H, G )
	local K, V, orderGDivOrderH, orderKIntHDivK, dimVKIntH, KIntH, conjecturedDim, obtainedDim, check;
	realIrr( H );
	orderGDivOrderH := Order( G )/Order( H );
	check := true;
	for K in AllSubgroups( G ) do
		KIntH := Intersection( K, H );
		orderKIntHDivK := Order( KIntH )/Order( K );
		for V in realIrreducibles do
			dimVKIntH := fixedPointDimensionRealModule( [[V,1]], KIntH, H );
			conjecturedDim := orderGDivOrderH*orderKIntHDivK*dimVKIntH;
			obtainedDim := fixedPointDimensionRealModuleInduced( [[V,1]], H, G, K );
			if conjecturedDim <> obtainedDim then
				#Print( "H = ", StructureDescription( H ), " = ", IdGroup( H ), " ",
				#	   "K = ", StructureDescription( K ), " = ", IdGroup( K ), "\n" );
				#Print( "Conjectured fixed point dimension: ", conjecturedDim, 
				#	   "\nObtained fixed point dimension:    ", obtainedDim, "\n" );
				check := false;
			fi;
		od;
	od;
	if check = true then
		Print( "H = ", StructureDescription( H ), " = ", IdGroup( H ), "\n" );
	fi;
end;

displayInductionFixedPointDimensionInfo := function( G )
	local H, i, subgroupsInductionInjective;
	subgroupsInductionInjective := [];
	for H in AllSubgroups( G ) do
		if isRealInductionHomomorphismInjective( H, G ) = true then
			Add( subgroupsInductionInjective, H );
		fi;
	od;
	i := 1;
	Print( "Number of subgroups with induction homomorphism injective: ", Size( subgroupsInductionInjective ), "\n" );
	for H in subgroupsInductionInjective do
		inducedDimensionsSubgroups( H, G );
	od;
end;

# Functions responsible for determining the groups with Smith sets equal
# to the reduced primary groups or contained in them.

# To load this file paste this line (change the directory if necessary),
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "smithSetPORed.g" ) );

# Requires "commonFunctions.g" to be called earlier. It can be called with the following command,
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "commonFunctions.g" ) );

# TO VERIFY
# Checks the sufficient conditions for a group G to have the Smith set contained
# in the reduced primary group. For details, see [1,UZUPELNIC] and [1,UZUPELNIC]
# Requires realIrr( G ) to be called earlier.
smithSetContainedInPORed := function( G )
	local g, W, check, H, h;
	for g in G do
		if IsPrimePowerInt( Order( g ) ) and ( Order( g ) mod 8 = 0 ) then
			check := 0;
			for W in realIrreducibles do
				if fixedPointDimensionRealModule( [[W,1]], GroupByGenerators( [g] ), G ) = 0 then
					check := check+1;
				fi;
			od;
			H := GroupByGenerators( [g^2] );
			for h in H do
				if IsConjugate( G, g*h, g ) = false and IsConjugate( G, g*h, g^(-1) ) = false then
					check := check+1;
					break;
				fi;
			od;
			if check = 2 then
				return false;
			fi;
		fi;
	od;
	return true;
end;

# TO VERIFY
# Checks the sufficient conditions for a group G to have the reduced primary group
# as the Smith set. For details, see [1,UZUPELNIC] and [1,UZUPELNIC]
# Requires realIrr( G ) to be called earlier.
hasPORedAsSmithSet := function( G )
	local g, W, check, H, h;
	if IsPerfect( G ) = false then
		return false;
	fi;
	return smithSetContainedInPORed( G );
end;

# REFERENCES
# [1] K. M. Pawałowski, The Smith equivalence problem and the Laitinen conjecture,
#	  in: Handbook of group actions. Vol III, volume 40. of Adv. Lect. Math. (ALM),
#	  Int. Press, Somerville, MA, 2018, pp. 485-537.

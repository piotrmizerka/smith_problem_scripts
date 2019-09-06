# To load this file separately, paste this line (change the directory if necessary),
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "constraintsSave.g" ) );

# Requires "commonFunctions.g", "genPORed.g", "latticePoints.sage", "realModulesGivenDimension.g"
# and "sufficientConditions.g" to be called earlier in that order.
# They can be called with the following commands,
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "commonFunctions.g" ) );
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "genPORed.g" ) );
# In terminal: sage --> load( "latticePoints.sage" )
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "realModulesGivenDimension.g" ) );
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "sufficientConditions.g" ) );

int2String := function( n )
  local result;
  result := "";
  while n > 0 do
    result := Concatenation( [CharInt( (n mod 10)+48 )], result );
    n := Int( n/10 );
  od;
  if result = "" then
    return "0";
  fi;
  return result;
end;

# Requires to be called earlier: realIrr( G ), saveGeneratorsPORed( G, "./group_data/gen_po_red.sage" ),
# determineConstraints( G ) and Read( "./group_data/lattice_points_red_po.g" ) (in GAP),
# saveLatticePoints( generatorVerticesFromFile, realIrreduciblesDimensionsFromFile, dimensionUpperBound ) (in SAGE).
constraintsSave := function( savePath, G, constraintsMinDimensions, minDimensionWToConsiderId )
  local latticePoint, i, currentSavePath;
  i := 1;
  for latticePoint in latticePointsRedPO do
    currentSavePath := savePath;
    currentSavePath := Concatenation( currentSavePath, "constraintId_" );
    currentSavePath := Concatenation( currentSavePath, int2String( i ) );
    currentSavePath := Concatenation( currentSavePath, "_minDimensionWToConsiderId_" );
    currentSavePath := Concatenation( currentSavePath, int2String( minDimensionWToConsiderId ) );
    currentSavePath := Concatenation( currentSavePath, ".sage" );
    saveConstraintsAsSAGEFile( currentSavePath, latticePoint[1], G, constraintsMinDimensions[i] );
    i := i+1;
  od;
end;

# To load this file separately, paste this line (change the directory if necessary),
# Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "specialOliver.g" ) );


isOliver := function( G )
  local N, H, idx;
  for H in NormalSubgroups( G ) do
    idx := IndexNC( G, H );
    if IsPrimePowerInt( idx ) or IsOne( idx ) then
      for N in NormalSubgroups( H ) do
        if IsPGroup( N ) and IsCyclic( H/N ) then
          return false;
        fi;
      od;
    fi;
  od;
  return true;
end;

hasNoNppOddOrderCyclicQuotient := function( G )
  local quotient, quotientOrder, N;
  for N in NormalSubgroups( G ) do
    quotient := G/N;
    quotientOrder := Order( quotient );
    if quotientOrder <> 1 and IsCyclic( quotient ) and (IsPrimePowerInt( quotientOrder )=false)
       and (IsEvenInt( quotientOrder ) = false) then
      Display( StructureDescription( quotient ) );
      return false;
    fi;
  od;
  return true;
end;

# TODO
satisfiesGnilCondition := function( G )
  return false;
end;

isSpecialOliver := function( G )
  return (isOliver( G ) and IsEvenInt( Order( G ) ) and hasNoNppOddOrderCyclicQuotient( G ) and
          IsSolvable( G ) and (satisfiesGnilCondition( G ) = false));
end;

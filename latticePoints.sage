# cd Desktop - > load( "latticePoints.sage" )

from sage.geometry.integral_points import simplex_points

# Step 1. function.
# Lattice points in a bounded convex neighborhood of the origin defined 
# by 2^n simplices defined by "generatorVertices" and "scale" (here n is the number of "generatorVertices").
# This region is a closed ball in the taxicab metric.
def latticePoints( generatorVertices, scale ):
	result = []
	simplexTuples = Tuples( [-1,1], len( generatorVertices ) ).list()
	origin = [0]*len( generatorVertices[0] )
	for tuple in simplexTuples:
		simplexVertices = []
		i = 0
		for vertex in generatorVertices:
			simplexVertex = []
			for coord in vertex:
				simplexVertex.append( coord*tuple[i]*scale )
			simplexVertices.append( simplexVertex )
			i += 1
		simplexVertices.append( origin )
		simplex = Polyhedron( simplexVertices )
		latticePointsx = simplex_points( simplex.Vrepresentation() )
		for point in latticePointsx:
			check = False
			for coord in point:
				if coord <> 0:
					check = True
					break
			if check == True:
				pointx = point
				pointx.set_immutable()
				result.append( pointx )
		resultx = set( result )
 		resultxx = list( resultx )
	return resultxx

# Step 2. function. 
# Computing the Eulcidean distance from the hyperplane spanned by "vertices" to the origin.
def distanceFromOriginHyperplane( vertices ):
	vectors = []
	for i in range( 1, len( vertices ) ):
		vectors.append( vertices[i]-vertices[0] )
	bas = list( (matrix( RDF,vectors ).gram_schmidt())[0] )
	return (vertices[0]-sum([(vertices[0]*x)*x for x in bas])).norm()

# Step 2. function. 
# Computing the Eulcidean distance from the boundary of the convex hull spanned by "generatorVertices"
def distanceFromOrigin( generatorVertices ):
	simplexTuples = Tuples( [-1,1], len( generatorVertices ) ).list()
	result = 1000000000
	for tuple in simplexTuples:
		simplexVertices = []
		i = 0
		for vertex in generatorVertices:
			simplexVertex = []
			for coord in vertex:
				simplexVertex.append( coord*tuple[i] )
			simplexVertices.append( vector( RDF,simplexVertex ) )
			i += 1
		dist = distanceFromOriginHyperplane( simplexVertices )
		if dist < result:
			result = dist
	return result

# Step 3. function.
# Finding the bounded region with at least one lattice point. 
# We look for it over all $\Delta(\lambda)$ for $\lambda=1,...,"maxScale"$
# Returns parameters "generatorVertices" and "i" for the bounded region 
# with lattice point if it was found together with the minimal lattice point found.
# This is why we need "irreduciblesDimensions" as parameter.
def boundedRegionLatticePoint( generatorVertices, maxScale, irreduciblesDimensions ):
	for i in range( 1, maxScale+1 ):
		latticePointsx = latticePoints( generatorVertices, i )
		if len( latticePointsx ) > 0:
			minIndex = 0
			minDoubledDimension = 1000000000
			for j in range( 0, len( latticePointsx ) ):
				doubledDimension = 0
				for k in range( 0, len( latticePointsx[j] ) ):
					doubledDimension += abs( latticePointsx[j][k] )*irreduciblesDimensions[k]
				if minDoubledDimension > doubledDimension:
					minDoubledDimension = doubledDimension
					minIndex = j
			return [generatorVertices,i,latticePointsx[minIndex]]

# Step 4. function.
# Enlarge the bounded region with lattice point so that it containes the minimal lattice point.
# Parameters: 
# - "generatorVertices" - generators of the reduced primary group
# - "scaleLatticePoint" - the scale for which a lattice point was found - output from "boundedRegionLatticePoint"
# - "minLatticePointFound" - minimal lattice point found by "boundedRegionLatticePoint"
# - "irreducibleDimensions" - dimensions of real irreducible representations
# Output:
# - the scale at which $\Delta(1)$ should be scaled to obtain the enlarged region
def enlargedRegion( generatorVertices, scaleLatticePoint, minLatticePointFound, irreduciblesDimensions ):
	latticePointFoundDoubledDimension = 0
	for i in range( 0, len( irreduciblesDimensions ) ):
		latticePointFoundDoubledDimension += abs( minLatticePointFound[i] )*irreduciblesDimensions[i]
	adjustedScale = int( latticePointFoundDoubledDimension/distanceFromOrigin( generatorVertices ) )+1
	return adjustedScale

# Step 5. function.
# Finds the minimal lattice point we look for. It calls all the necessary functions defined earlier.
# Parameters:
# - "generatorVertices" - generators of the reduced primary group
# - "irreduciblesDimensions" - dimensions of real irreducible representations
# - "maxScale" - the scale up to which look for the lattice point
# Output:
# - minimal lattice point U-V and the dimension of RG-modules U and V
def minimalLatticePoint( generatorVertices, irreduciblesDimensions, maxScale = 100 ):
	[generatorVerticesx,scalex,latticePointx] = boundedRegionLatticePoint( generatorVertices, maxScale, 
                                                                               irreduciblesDimensions )
	adjustedScale = enlargedRegion( generatorVertices, scalex, latticePointx, irreduciblesDimensions )
	latticePointsx = latticePoints( generatorVertices, adjustedScale )
	minNorm = 1000000000
	minIndex = 0
	for i in range( 0, len( latticePointsx ) ):
		norm = 0
		for j in range( 0, len( irreduciblesDimensions ) ):
			norm += abs( latticePointsx[i][j] )*irreduciblesDimensions[j]
		if norm < minNorm:
			minNorm = norm
			minIndex = i
	return [latticePointsx[minIndex],minNorm/2]
	#return adjustedScale

# Computes all the lattice points U-V for which dimU=dimV<="dimensionUpperBound".
# Parameters:
# - "generatorVertices" - generators of the reduced primary group
# - "irreduciblesDimensions" - dimensions of real irreducible representtions
# - "dimensionUpperBound" - upper bound for dimensions to look for
# Output:
# - lattice points together with dimensions they define in the format [dimension,lattice point]
def latticePointsUpToDimension( generatorVertices, irreduciblesDimensions, dimensionUpperBound ):
	scale = int( (2*dimensionUpperBound)/distanceFromOrigin( generatorVertices ) )+1
	latticePointsx = latticePoints( generatorVertices, scale )
	latticePointsxx = []
	for latticePoint in latticePointsx:
		norm = 0
		for i in range( 0, len( latticePoint ) ):
			norm += abs( latticePoint[i] )*irreduciblesDimensions[i]
		if norm <= 2*dimensionUpperBound:
			latticePointsxx.append( (norm/2,latticePoint) )
	latticePointsxxx = sorted( latticePointsxx, key = lambda y: y[0] )
	return latticePointsxx
    
# Saves to the gap file all the lattice points U-V for which dimU=dimV<="dimensionUpperBound".
# These points may be used later by an appropriate GAP-procedure to establish which of them satisfy 
# the conditions providing existince of an exotic Smith action.
# Parameters:
# - "generatorVertices" - generators of the reduced primary group
# - "irreduciblesDimensions" - dimensions of real irreducible representtions
# - "dimensionUpperBound" - upper bound for dimensions to look for
def saveLatticePoints( generatorVertices, irreduciblesDimensions, dimensionUpperBound ):
	latticePointsxxx = latticePointsUpToDimension( generatorVertices, irreduciblesDimensions, dimensionUpperBound )
	save = open( "smith_set_elements.g", "w" )
	save.write( "smithSetElements := [" )
	for i in range( len( latticePointsxxx ) ):
		if i > 0:
			save.write( "," )
		save.write( "[[" )
		for j in range( 0, len( latticePointsxxx[i][1] )-1 ):
			save.write( str( latticePointsxxx[i][1][j] )+"," )
		save.write( str( latticePointsxxx[i][1][len( latticePointsxxx[i][1] )-1] )+"],"+str( latticePointsxxx[i][0] )+"]" )
	save.write( "];" )
	save.close()
			
#lattice_points = latticePoints( [(1,0),(0,1)],3 )
#print( len( lattice_points ) )
#print( distanceFromOrigin( [(1,2,-10),(2,1,-10)] ) )
#print( distanceFromOriginSimplex( [vector( RDF,[-2,-1,10] ), vector( RDF,[1,2,-10] )] ) )
#print( minimalLatticePoint( [(1,-1,-2,2,0,0,0,0),(1,0,-2,0,4,-2,-2,0)],(4,4,3,3,4,8,5,12) ) ) # SL(2,5) - lower bound: 10
#print( minimalLatticePoint( [(1,-1,2,-1,-1,2,-4,2)],(6,8,6,12,12,7,8,16) ) ) # SL(2,7) - lower bound: 64
#print( minimalLatticePoint( [(2,0,0,-2,-1,1,0,0,0,0,0,0),(2,0,0,-2,-1,0,1,0,0,0,0,0),(-11,-5,-8,8,3,0,0,10,0,0,0,0),
#(-12,0,-6,6,6,0,0,0,-5,0,10,0),(-12,0,-6,6,6,0,0,0,0,-5,0,10)],(10,12,10,10,20,20,20,11,24,24,12,12) ) ) # SL(2,11) - lower bound: 40
#saveLatticePoints( [(1,-1,-2,2,0,0,0,0),(1,0,-2,0,4,-2,-2,0)],(4,4,3,3,4,8,5,12), 30 )
#saveLatticePoints( [(1,-1,2,-1,-1,2,-4,2)],(6,8,6,12,12,7,8,16), 90 )
#saveLatticePoints( [(2,0,0,-2,-1,1,0,0,0,0,0,0),(2,0,0,-2,-1,0,1,0,0,0,0,0),(-11,-5,-8,8,3,0,0,10,0,0,0,0),
#(-12,0,-6,6,6,0,0,0,-5,0,10,0),(-12,0,-6,6,6,0,0,0,0,-5,0,10)],(10,12,10,10,20,20,20,11,24,24,12,12), 40 )
	
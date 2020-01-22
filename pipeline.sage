# To load this file run sage and call "load( "pipeline.sage" )"

import prettytable
from prettytable import PrettyTable

# Global variables necessary for running the program
objectiveValuex = [0]*100
objectiveCoordinatesx = [0]*100
suitableSmithPairsIds = [False]*100

def saveConstraints( minDimensionWToConsiderId, constraintsMinDimensions, latticePointsRedPONumber ):
    constraintsMinDimensions = [0]*latticePointsRedPONumber
    if minDimensionWToConsiderId == 0:
        gap.eval( 'constraintsSave( "./group_data/constraints/", G, '+str( constraintsMinDimensions )+' , 0 );' )
    else:
        for i in range( latticePointsRedPONumber ):
            load( "./group_data/constraints/constraintId_"+str( i+1 )+"_minDimensionWToConsiderId_"+str( minDimensionWToConsiderId-1 )+".sage" )
            constraintsMinDimensions[i] = int( objectiveValue )+1
        gap.eval( 'constraintsSave( "./group_data/constraints/", G, '+str( constraintsMinDimensions )+', '+str( minDimensionWToConsiderId )+' );' )

def determineSuitableSmithPairsIds( latticePointsRedPONumber, amountOfOptimizedPointsToCheck, realIrreduciblesNumber ):
    for i in range( 1, latticePointsRedPONumber+1 ):
        suitableSmithPairsIds[i] = [False]*amountOfOptimizedPointsToCheck
        for j in range( amountOfOptimizedPointsToCheck ):
            gapEvalStr = 'U := []; V := []; W := [];\n'
            gapEvalStr += 'for i in [1..Size( realIrreducibles )] do\n'
            gapEvalStr += '\tif latticePointsRedPO['+str( i )+'][1][i] > 0 then\n'
            gapEvalStr += '\t\tAdd( U, [realIrreducibles[i],'+'latticePointsRedPO['+str( i )+'][1][i]'+'] );\n'
            gapEvalStr += '\telif latticePointsRedPO['+str( i )+'][1][i] < 0 then\n'
            gapEvalStr += '\t\tAdd( V, [realIrreducibles[i],'+'-latticePointsRedPO['+str( i )+'][1][i]'+'] );\n'
            gapEvalStr += '\tfi;\n'
            gapEvalStr += 'od;\n'
            for k in range( realIrreduciblesNumber ):
                if objectiveCoordinatesx[i][j][k][1] > 0:
                    gapEvalStr += ("Add( W, [realIrreducibles["+str( objectiveCoordinatesx[i][j][k][0]+1 )+"],"+str( int( objectiveCoordinatesx[i][j][k][1] ) )+"] );")
            gap.eval( gapEvalStr )
            isPOrientedPair = gap.eval( 'verifyPOrientabilityModule( U, W, G ) and verifyPOrientabilityModule( V, W, G )' )
            if isPOrientedPair:
                suitableSmithPairsIds[i][j] = True

def saveSuitableSmithPairs( latticePointsRedPONumber, amountOfOptimizedPointsToCheck, realIrreduciblesNumber ):
    latticePointsRedPOx = gap.eval( 'latticePointsRedPO' )
    latticePointsRedPOx = "latticePointsRedPOAsSageVariable = "+latticePointsRedPOx
    exec( latticePointsRedPOx )
    realIrreduciblesx = gap.eval( 'realIrreducibles' )
    realIrreduciblesx = realIrreduciblesx.replace( '^', '**' )
    realIrreduciblesx = "realIrreduciblesAsSageVariable = "+realIrreduciblesx
    exec( realIrreduciblesx )
    groupAsString = gap.eval( 'G' )
    save = open( "./group_data/"+groupAsString+"_Smith_exoticism.txt", "w" )
    save.write( '# Complex character table\n\n' )
    save.write( gap.eval('Display(CharacterTable(G))')+'\n' )
    save.write( '\n\n# Real irreducible characters\n\n' )
    realCharacterTable = PrettyTable()
    realCharacterTable.hrules = prettytable.FRAME
    temp = gap.eval( 'for cl in ConjugacyClasses(G) do Print(Order(Representative(cl))," ");od;' ).split()
    d = {}
    for i in range( len( temp ) ):
        d[temp[i]] = 0
    tableRow = []
    tableRow.append( 'conjugacy_class' )
    for i in range( len( temp ) ):
        tableRow.append( str( temp[i] )+chr( d[temp[i]]+97 ) )
        d[temp[i]] += 1
    realCharacterTable.field_names = tableRow
    tableRow = 'conjugacy_class_size '+gap.eval( 'for cl in ConjugacyClasses(G) do Print(Size(cl)," ");od;' )
    tableRow = tableRow.split()
    realCharacterTable.add_row( tableRow )
    for i in range( realIrreduciblesNumber ):
        realCharacterTable.add_row( ['V.'+str( i+1 )]+realIrreduciblesAsSageVariable[i] )
    save.write( str( realCharacterTable )+"\n" )
    save.write( '\n\n# Fixed point dimension table for real irreducible representation\n\n' )
    exec( 'tableFixedPointDimensionSage = '+gap.eval( "tableFixedPointDimension( G )" ) )
    tableFixedPointDimension = PrettyTable()
    tableFixedPointDimension.hrules = prettytable.ALL
    tableFixedPointDimension.field_names = ["group_class_id"]+tableFixedPointDimensionSage[0]
    for i in range( realIrreduciblesNumber ):
        tableFixedPointDimension.add_row( ["V."+str( i )]+tableFixedPointDimensionSage[i+1] )
    save.write( str( tableFixedPointDimension )+"\n" )
    save.write( '\n\n# Suitable Smith pairs\n\n' )
    suitableSmithPairsTable = PrettyTable()
    suitableSmithPairsTable.hrules = prettytable.ALL
    suitableSmithPairsTable.field_names = ['U','V','dimension']
    load( './group_data/gen_po_red.sage' )
    dimensionsGcd = 0
    for i in range( latticePointsRedPONumber ):
        for j in range( amountOfOptimizedPointsToCheck ):
            U = ''
            V = ''
            dimension = 0
            if suitableSmithPairsIds[i+1][j]:
                checkU = 0
                checkV = 0
                for k in range( realIrreduciblesNumber ):
                    coefficientU = 0
                    coefficientV = 0
                    if latticePointsRedPOAsSageVariable[i][0][k] > 0:
                        coefficientU = latticePointsRedPOAsSageVariable[i][0][k]+int( objectiveCoordinatesx[i+1][j][k][1] )
                        coefficientV = int( objectiveCoordinatesx[i+1][j][k][1] )
                        checkU += 1
                    elif latticePointsRedPOAsSageVariable[i][0][k] < 0:
                        coefficientU = int( objectiveCoordinatesx[i+1][j][k][1] )
                        coefficientV = -latticePointsRedPOAsSageVariable[i][0][k]+int( objectiveCoordinatesx[i+1][j][k][1] )
                        checkV += 1
                    else:
                        coefficientU = int( objectiveCoordinatesx[i+1][j][k][1] )
                        coefficientV = int( objectiveCoordinatesx[i+1][j][k][1] )
                    if coefficientU > 0:
                        dimension += coefficientU*realIrreduciblesDimensionsFromFile[k]
                        if checkU > 1:
                            U += '+'
                        if coefficientU > 1:
                            U += (str( coefficientU )+'V.'+str( k+1 ))
                        else:
                            U += ('V.'+str( k+1 ))
                    if coefficientV > 0:
                        if checkV > 1:
                            V += '+'
                        if coefficientV > 1:
                            V += (str( coefficientV )+'V.'+str( k+1 ))
                        else:
                            V += ('V.'+str( k+1 ))
            suitableSmithPairsTable.add_row( [U,V,dimension] )
            if dimensionsGcd == 0:
                dimensionsGcd = dimension
            else:
                dimensionsGcd = gcd( dimensionsGcd, dimension )
    save.write( str( suitableSmithPairsTable )+"\n" )
    save.write( "\n\n# Dimensions gretest common divisor\n\n"+str( dimensionsGcd ) )
    save.close()

def pipeline( groupAsString, amountOfOptimizedPointsToCheck ):
    gap.eval( 'Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "commonFunctions.g" ) );' )
    gap.eval( 'Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "genPORed.g" ) );' )
    gap.eval( 'G := '+groupAsString+';' )
    gap.eval( 'realIrr( G );' )
    gap.eval( 'saveGeneratorsPORed( G, "./group_data/gen_po_red.sage" );' )
    load( './group_data/gen_po_red.sage' )
    load( 'latticePoints.sage' )
    saveLatticePoints( genaratorVerticesFromFile, realIrreduciblesDimensionsFromFile, 30 )
    gap.eval( 'Read( "./group_data/lattice_points_red_po.g" );' )
    gap.eval( 'Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "realModulesGivenDimension.g" ) );' )
    gap.eval( 'Read( Filename( [DirectoryCurrent(), DirectoryDesktop()], "sufficientConditions.g" ) );' )
    gap.eval( 'determineConstraints( G );' )
    gap.eval( 'Read( "constraintsSave.g" );' )
    latticePointsRedPONumber = int( gap.eval( 'Size( latticePointsRedPO )' ) )
    realIrreduciblesNumber = int( gap.eval( 'Size( realIrreducibles )' ) )
    constraintsMinDimensions = [0]*latticePointsRedPONumber
    for j in range( amountOfOptimizedPointsToCheck ):
        saveConstraints( j, constraintsMinDimensions, latticePointsRedPONumber )
    for i in range( 100 ):
        objectiveValuex[i] = [0]*100
        objectiveCoordinatesx[i] = [[]]*100
        suitableSmithPairsIds[i] = [False]*100
    for i in range( 1, latticePointsRedPONumber+1 ):
        for j in range( amountOfOptimizedPointsToCheck ):
            load( "./group_data/constraints/constraintId_"+str( i )+"_minDimensionWToConsiderId_"+str( j )+".sage" )
            objectiveValuex[i][j] = objectiveValue
            objectiveCoordinatesx[i][j] = objectiveCoordinates
    determineSuitableSmithPairsIds( latticePointsRedPONumber, amountOfOptimizedPointsToCheck, realIrreduciblesNumber )
    saveSuitableSmithPairs( latticePointsRedPONumber, amountOfOptimizedPointsToCheck, realIrreduciblesNumber )

# SAMPLE USAGE
#pipeline( "SL(2,5)", 3 )

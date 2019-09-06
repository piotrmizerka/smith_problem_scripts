equationSystem = MixedIntegerLinearProgram( maximization=False, solver = "GLPK\
" )
x = equationSystem.new_variable( integer = True, nonnegative = True )


#Constraints from gap hypothesis
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(-3)*x[2]+(-3)*x[3]+(-4)*x[4]+(8)*x[5]+(-5)*x[6]+(12)*x[7] >= 2 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(-3)*x[2]+(-3)*x[3]+(-4)*x[4]+(8)*x[5]+(-5)*x[6]+(12)*x[7] >= 2 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(3)*x[6]+(4)*x[7] >= -6 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(3)*x[6]+(4)*x[7] >= -6 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(8)*x[5]+(-1)*x[6]+(12)*x[7] >= -6 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(8)*x[5]+(-1)*x[6]+(12)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(1)*x[2]+(1)*x[3]+(4)*x[4]+(8)*x[5]+(3)*x[6]+(4)*x[7] >= -6 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(1)*x[2]+(1)*x[3]+(4)*x[4]+(8)*x[5]+(3)*x[6]+(4)*x[7] >= -6 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(8)*x[5]+(3)*x[6]+(12)*x[7] >= -6 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(8)*x[5]+(3)*x[6]+(12)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(-1)*x[2]+(-1)*x[3]+(-2)*x[4]+(4)*x[5]+(-1)*x[6]+(4)*x[7] >= 2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(-1)*x[2]+(-1)*x[3]+(-2)*x[4]+(4)*x[5]+(-1)*x[6]+(4)*x[7] >= 2 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(8)*x[5]+(1)*x[6]+(12)*x[7] >= -10 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(8)*x[5]+(1)*x[6]+(12)*x[7] >= -10 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(1)*x[2]+(1)*x[3]+(4)*x[4]+(8)*x[5]+(3)*x[6]+(12)*x[7] >= -6 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(1)*x[2]+(1)*x[3]+(4)*x[4]+(8)*x[5]+(3)*x[6]+(12)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(4)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(4)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(-1)*x[2]+(-1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(4)*x[7] >= 2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(-1)*x[2]+(-1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(4)*x[7] >= 2 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(8)*x[5]+(3)*x[6]+(12)*x[7] >= -10 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(8)*x[5]+(3)*x[6]+(12)*x[7] >= -10 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(4)*x[5]+(-1)*x[6]+(4)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(4)*x[5]+(-1)*x[6]+(4)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(3)*x[6]+(12)*x[7] >= -10 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(3)*x[6]+(12)*x[7] >= -10 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(4)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(-1)*x[6]+(4)*x[7] >= -2 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -10 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -10 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(2)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(-1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(-1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 0 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -10 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -10 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -6 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(4)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(4)*x[7] >= -2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 0 )

#Constraints from condition 2
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -5 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= 3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= 3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= 3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= 3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(4)*x[7] >= 3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(4)*x[7] >= 3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 5 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -8 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -8 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -4 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -4 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(4)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(4)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 0 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(0)*x[6]+(0)*x[7] >= 2 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(0)*x[6]+(0)*x[7] >= 2 )

#Constraints from condition 3
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(8)*x[5]+(0)*x[6]+(12)*x[7] >= -3 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(8)*x[5]+(0)*x[6]+(12)*x[7] >= -3 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(2)*x[2]+(2)*x[3]+(2)*x[4]+(4)*x[5]+(4)*x[6]+(8)*x[7] >= -7 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(2)*x[2]+(2)*x[3]+(2)*x[4]+(4)*x[5]+(4)*x[6]+(8)*x[7] >= -7 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(2)*x[2]+(2)*x[3]+(2)*x[4]+(8)*x[5]+(2)*x[6]+(12)*x[7] >= -7 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(2)*x[2]+(2)*x[3]+(2)*x[4]+(8)*x[5]+(2)*x[6]+(12)*x[7] >= -7 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(2)*x[2]+(2)*x[3]+(2)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= -3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(2)*x[2]+(2)*x[3]+(2)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= -3 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(2)*x[2]+(2)*x[3]+(4)*x[4]+(8)*x[5]+(4)*x[6]+(8)*x[7] >= -7 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(2)*x[2]+(2)*x[3]+(4)*x[4]+(8)*x[5]+(4)*x[6]+(8)*x[7] >= -7 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(2)*x[2]+(2)*x[3]+(2)*x[4]+(8)*x[5]+(4)*x[6]+(12)*x[7] >= -7 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(2)*x[2]+(2)*x[3]+(2)*x[4]+(8)*x[5]+(4)*x[6]+(12)*x[7] >= -7 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(2)*x[2]+(2)*x[3]+(2)*x[4]+(0)*x[5]+(4)*x[6]+(0)*x[7] >= -3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(2)*x[2]+(2)*x[3]+(2)*x[4]+(0)*x[5]+(4)*x[6]+(0)*x[7] >= -3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(4)*x[5]+(0)*x[6]+(4)*x[7] >= 1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(4)*x[5]+(0)*x[6]+(4)*x[7] >= 1 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(8)*x[5]+(3)*x[6]+(12)*x[7] >= -9 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(8)*x[5]+(3)*x[6]+(12)*x[7] >= -9 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(2)*x[2]+(2)*x[3]+(4)*x[4]+(8)*x[5]+(4)*x[6]+(12)*x[7] >= -7 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(2)*x[2]+(2)*x[3]+(4)*x[4]+(8)*x[5]+(4)*x[6]+(12)*x[7] >= -7 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(2)*x[2]+(2)*x[3]+(4)*x[4]+(0)*x[5]+(4)*x[6]+(0)*x[7] >= -3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(2)*x[2]+(2)*x[3]+(4)*x[4]+(0)*x[5]+(4)*x[6]+(0)*x[7] >= -3 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(0)*x[5]+(0)*x[6]+(4)*x[7] >= 1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(0)*x[5]+(0)*x[6]+(4)*x[7] >= 1 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(8)*x[5]+(4)*x[6]+(12)*x[7] >= -9 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(8)*x[5]+(4)*x[6]+(12)*x[7] >= -9 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(0)*x[5]+(4)*x[6]+(0)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(0)*x[5]+(4)*x[6]+(0)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(4)*x[5]+(0)*x[6]+(4)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(4)*x[5]+(0)*x[6]+(4)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(0)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(0)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(4)*x[6]+(12)*x[7] >= -9 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(4)*x[6]+(12)*x[7] >= -9 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(4)*x[6]+(0)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(4)*x[6]+(0)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(0)*x[6]+(4)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(0)*x[6]+(4)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(0)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(0)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -9 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -9 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(3)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 1 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -9 )
equationSystem.add_constraint( (4)*x[0]+(4)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(8)*x[5]+(5)*x[6]+(12)*x[7] >= -9 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(3)*x[2]+(3)*x[3]+(4)*x[4]+(0)*x[5]+(5)*x[6]+(0)*x[7] >= -5 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(4)*x[5]+(1)*x[6]+(4)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(3)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(4)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(4)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(2)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(2)*x[6]+(0)*x[7] >= 1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(1)*x[2]+(1)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= -1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(0)*x[4]+(0)*x[5]+(1)*x[6]+(0)*x[7] >= 1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(0)*x[6]+(0)*x[7] >= 1 )
equationSystem.add_constraint( (0)*x[0]+(0)*x[1]+(0)*x[2]+(0)*x[3]+(1)*x[4]+(0)*x[5]+(0)*x[6]+(0)*x[7] >= 1 )

#Constraints from condition 4

#Constraint for minimal dimension of W
equationSystem.add_constraint( 4*x[0]+4*x[1]+3*x[2]+3*x[3]+4*x[4]+8*x[5]+5*x[6]+12*x[7] >= 43 )

equationSystem.set_objective( 4*x[0]+4*x[1]+3*x[2]+3*x[3]+4*x[4]+8*x[5]+5*x[6]+12*x[7] )

# Results of the above integer linear porgramming problem
objectiveValue = equationSystem.solve()
objectiveCoordinates = sorted( equationSystem.get_values( x ).items() )

#Uncomment the following lines for more desciptive output,
#equationSystem.show()

#print( 'Objective Value: {}'.format( equationSystem.solve() ) )
#for i, v in sorted( equationSystem.get_values( x ).items()):
#	print( 'w_%s = %s' % (i, int( round( v ) )) )
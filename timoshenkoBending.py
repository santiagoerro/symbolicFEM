import sympy as sym


l = sym.Symbol('l', positive = True)
shearCorrectionFactorVertical = sym.Symbol('s_v', positive = True)
shearCorrectionFactorHorizontal = sym.Symbol('s_h', positive = True)
a = sym.Symbol('a', positive = True)
EA = sym.Symbol('EA')
EIVertical = sym.Symbol('EI_v')
EIHorizontal = sym.Symbol('EI_h')
linearDensity = sym.Symbol('lambda')
zCenterOfMass = sym.Symbol('z_m')
zNeutralAxis = sym.Symbol('z_n')
zCenterOfTwist = sym.Symbol('z_t')
verticalMassInertiaMoment = sym.Symbol('I_yy')
horizontalMassInertiaMoment = sym.Symbol('I_zz')
rollMassInertiaMoment = sym.Symbol('I_xx')

# C11w0 = sym.Symbol('C11w0')
# C12w0 = sym.Symbol('C12w0')
# C13w0 = sym.Symbol('C13w0')
# C14w0 = sym.Symbol('C14w0')
# C11w1 = sym.Symbol('C11w1')
# C12w1 = sym.Symbol('C12w1')
# C13w1 = sym.Symbol('C13w1')
# C14w1 = sym.Symbol('C14w1')

x = sym.Symbol('x')

printType = 'pretty'



def SimplifyPrint(expression, factor = None, series = False, symbol = sym.Symbol('x'), x0 = 0, order = 10):
    if printType == 'latex':
        Printer = lambda x: print(sym.latex(x))
    elif printType == 'python':
        Printer = lambda x: print(sym.python(x))
    else:
        Printer = sym.pprint

    if factor is None:
        simplifiedExpression = sym.simplify(expression)
        Printer(simplifiedExpression)
        if series:
            Printer(sym.series(simplifiedExpression, symbol, x0, order))
    else:
        simplifiedFactor = sym.simplify(factor)
        simplifiedFactoredExpression = sym.simplify(expression / factor)
        Printer(simplifiedFactor)
        print('*')
        Printer(simplifiedFactoredExpression)
        if series:
            Printer(simplifiedFactor)
            print('*')
            Printer(sym.series(simplifiedFactoredExpression, symbol, x0, order))



axialDisplacementDisplacementBasis: sym.Matrix = sym.Matrix([1 - x/l, x/l]).transpose()

axialForceDisplacementBasis: sym.Matrix = EA * sym.diff(axialDisplacementDisplacementBasis, x)

axialStiffnessMatrix: sym.Matrix = sym.zeros(2, 2)
axialStiffnessMatrix[0, :] = -axialForceDisplacementBasis.subs(x, 0)
axialStiffnessMatrix[1, :] = axialForceDisplacementBasis.subs(x, l)


verticalDeflectionMonomialBasis: sym.Matrix = sym.Matrix([1, x, x**2, x**3]).transpose()
pitchRotationMonomialBasis: sym.Matrix = -sym.diff(verticalDeflectionMonomialBasis, x) - shearCorrectionFactorVertical * l**2 / 12 * sym.diff(verticalDeflectionMonomialBasis, x, 3)

horizontalDeflectionMonomialBasis: sym.Matrix = sym.Matrix([1, x, x**2, x**3]).transpose()
yawRotationMonomialBasis: sym.Matrix = sym.diff(verticalDeflectionMonomialBasis, x) + shearCorrectionFactorHorizontal * l**2 / 12 * sym.diff(verticalDeflectionMonomialBasis, x, 3)

monomialBasisCoefsToVerticalDisplacementsMatrix: sym.Matrix = sym.zeros(4, 4)
monomialBasisCoefsToVerticalDisplacementsMatrix[0, :] = verticalDeflectionMonomialBasis.subs(x, 0)
monomialBasisCoefsToVerticalDisplacementsMatrix[1, :] = pitchRotationMonomialBasis.subs(x, 0)
monomialBasisCoefsToVerticalDisplacementsMatrix[2, :] = verticalDeflectionMonomialBasis.subs(x, l)
monomialBasisCoefsToVerticalDisplacementsMatrix[3, :] = pitchRotationMonomialBasis.subs(x, l)

monomialBasisCoefsToHorizontalDisplacementsMatrix: sym.Matrix = sym.zeros(4, 4)
monomialBasisCoefsToHorizontalDisplacementsMatrix[0, :] = horizontalDeflectionMonomialBasis.subs(x, 0)
monomialBasisCoefsToHorizontalDisplacementsMatrix[1, :] = yawRotationMonomialBasis.subs(x, 0)
monomialBasisCoefsToHorizontalDisplacementsMatrix[2, :] = horizontalDeflectionMonomialBasis.subs(x, l)
monomialBasisCoefsToHorizontalDisplacementsMatrix[3, :] = yawRotationMonomialBasis.subs(x, l)

verticalDisplacementsToMonomialBasisCoefsMatrix: sym.Matrix = monomialBasisCoefsToVerticalDisplacementsMatrix.inv()
verticalDeflectionDisplacementBasis: sym.Matrix = verticalDeflectionMonomialBasis * verticalDisplacementsToMonomialBasisCoefsMatrix
pitchRotationDisplacementBasis: sym.Matrix = pitchRotationMonomialBasis * verticalDisplacementsToMonomialBasisCoefsMatrix

horizontalDisplacementsToMonomialBasisCoefsMatrix: sym.Matrix = monomialBasisCoefsToHorizontalDisplacementsMatrix.inv()
horizontalDeflectionDisplacementBasis: sym.Matrix = horizontalDeflectionMonomialBasis * horizontalDisplacementsToMonomialBasisCoefsMatrix
yawRotationDisplacementBasis: sym.Matrix = yawRotationMonomialBasis * horizontalDisplacementsToMonomialBasisCoefsMatrix

verticalBendingMomentDisplacementBasis: sym.Matrix = EIVertical * sym.diff(pitchRotationDisplacementBasis, x)
verticalShearForceDisplacementBasis: sym.Matrix = 12 * EIVertical / (shearCorrectionFactorVertical * l**2) * (sym.diff(verticalDeflectionDisplacementBasis, x) + pitchRotationDisplacementBasis)

horizontalBendingMomentDisplacementBasis: sym.Matrix = EIHorizontal * sym.diff(yawRotationDisplacementBasis, x)
horizontalShearForceDisplacementBasis: sym.Matrix = 12 * EIHorizontal / (shearCorrectionFactorHorizontal * l**2) * (sym.diff(horizontalDeflectionDisplacementBasis, x) - yawRotationDisplacementBasis)


verticalStiffnessMatrix: sym.Matrix = sym.zeros(4, 4)
verticalStiffnessMatrix[0, :] = -verticalShearForceDisplacementBasis.subs(x, 0)
verticalStiffnessMatrix[1, :] = -verticalBendingMomentDisplacementBasis.subs(x, 0)
verticalStiffnessMatrix[2, :] = verticalShearForceDisplacementBasis.subs(x, l)
verticalStiffnessMatrix[3, :] = verticalBendingMomentDisplacementBasis.subs(x, l)

horizontalStiffnessMatrix: sym.Matrix = sym.zeros(4, 4)
horizontalStiffnessMatrix[0, :] = -horizontalShearForceDisplacementBasis.subs(x, 0)
horizontalStiffnessMatrix[1, :] = -horizontalBendingMomentDisplacementBasis.subs(x, 0)
horizontalStiffnessMatrix[2, :] = horizontalShearForceDisplacementBasis.subs(x, l)
horizontalStiffnessMatrix[3, :] = horizontalBendingMomentDisplacementBasis.subs(x, l)


axialMassMatrix: sym.Matrix = linearDensity * sym.integrate(axialDisplacementDisplacementBasis.transpose() * axialDisplacementDisplacementBasis, (x, 0, l))
verticalBendingHeaveContributionMassMatrix: sym.Matrix = linearDensity * sym.integrate(verticalDeflectionDisplacementBasis.transpose() * verticalDeflectionDisplacementBasis, (x, 0, l))
verticalBendingPitchContributionMassMatrix: sym.Matrix = (verticalMassInertiaMoment + linearDensity * (zCenterOfMass - zNeutralAxis)**2) * sym.integrate(pitchRotationDisplacementBasis.transpose() * pitchRotationDisplacementBasis, (x, 0, l))
axialVerticalBendingCouplingMassMatrix: sym.Matrix = linearDensity * (zCenterOfMass - zNeutralAxis) * sym.integrate(axialDisplacementDisplacementBasis.transpose() * pitchRotationDisplacementBasis, (x, 0, l))
horizontalBendingSwayContributionMassMatrix: sym.Matrix = linearDensity * sym.integrate(horizontalDeflectionDisplacementBasis.transpose() * horizontalDeflectionDisplacementBasis, (x, 0, l))
horizontalBendingYawContributionMassMatrix: sym.Matrix = horizontalMassInertiaMoment * sym.integrate(yawRotationDisplacementBasis.transpose() * yawRotationDisplacementBasis, (x, 0, l))


twistingBoundedBasis: sym.Matrix = sym.Matrix([1, x/l, sym.exp(a * (x/l - 1)), sym.exp(-a * x/l)]).transpose()

boundedBasisCoefsToConditionsMatrix: sym.Matrix = sym.zeros(4, 4)
boundedBasisCoefsToConditionsMatrix[0, :] = twistingBoundedBasis.subs(x, 0)
boundedBasisCoefsToConditionsMatrix[1, :] = twistingBoundedBasis.subs(x, l)
boundedBasisCoefsToConditionsMatrix[2, :] = sym.diff(twistingBoundedBasis, x).subs(x, 0) * l
boundedBasisCoefsToConditionsMatrix[3, :] = sym.diff(twistingBoundedBasis, x).subs(x, l) * l
stableBasisToBoundedBasisCoefsMatrix: sym.Matrix = sym.zeros(4,4)
stableBasisToBoundedBasisCoefsMatrix[:, 0] = sym.Matrix([1, -1, 0, 0])
stableBasisToBoundedBasisCoefsMatrix[:, 1] = sym.simplify(boundedBasisCoefsToConditionsMatrix.LDLsolve(sym.Matrix([0, 0, 1, 0])))
stableBasisToBoundedBasisCoefsMatrix[:, 2] = sym.Matrix([0, 1, 0, 0])
stableBasisToBoundedBasisCoefsMatrix[:, 3] = sym.simplify(boundedBasisCoefsToConditionsMatrix.LDLsolve(sym.Matrix([0, 0, 0, 1])))

# stableBasisToBoundedBasisCoefsMatrix: sym.Matrix = sym.zeros(4,4)
# stableBasisToBoundedBasisCoefsMatrix[:, 0] = sym.Matrix([1, -1, 0, 0])
# stableBasisToBoundedBasisCoefsMatrix[:, 1] = sym.Matrix([C11w0, C12w0, C13w0, C14w0])
# stableBasisToBoundedBasisCoefsMatrix[:, 2] = sym.Matrix([0, 1, 0, 0])
# stableBasisToBoundedBasisCoefsMatrix[:, 3] = sym.Matrix([C11w1, C12w1, C13w1, C14w1])

massMatrixBoundedBasis = (rollMassInertiaMoment + linearDensity * (zCenterOfMass - zCenterOfTwist)**2) * sym.integrate(twistingBoundedBasis.transpose() * twistingBoundedBasis, (x, 0, l))
massMatrixStableBasis = sym.simplify(stableBasisToBoundedBasisCoefsMatrix.transpose() * massMatrixBoundedBasis * stableBasisToBoundedBasisCoefsMatrix)

horizontalBendingDispBasisTorsionBoundedBasisCouplingMassMatrix: sym.Matrix = -linearDensity * (zCenterOfMass - zCenterOfTwist) * sym.integrate(horizontalDeflectionDisplacementBasis.transpose() * twistingBoundedBasis, (x, 0, l))
horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix: sym.Matrix = sym.simplify(horizontalBendingDispBasisTorsionBoundedBasisCouplingMassMatrix * stableBasisToBoundedBasisCoefsMatrix)


print()
print('Axial force:')
SimplifyPrint(axialForceDisplacementBasis, factor = EA/l)
print()
print('Axial stiffness matrix:')
SimplifyPrint(axialStiffnessMatrix, factor = EA/l)
print()
print('Vertical bending moment:')
SimplifyPrint(verticalBendingMomentDisplacementBasis, factor = EIVertical / l**2 / (1 + shearCorrectionFactorVertical))
print()
print('Horizontal bending moment:')
SimplifyPrint(horizontalBendingMomentDisplacementBasis, factor = EIHorizontal / l**2 / (1 + shearCorrectionFactorHorizontal))
print()
print('Vertical shear force:')
SimplifyPrint(verticalShearForceDisplacementBasis, factor = EIVertical / l**3 / (1 + shearCorrectionFactorVertical))
print()
print('Horizontal shear force:')
SimplifyPrint(horizontalShearForceDisplacementBasis, factor = EIHorizontal / l**3 / (1 + shearCorrectionFactorHorizontal))
print()
print('Vertical bending stiffness matrix:')
SimplifyPrint(verticalStiffnessMatrix, factor = EIVertical / l**2 / (1 + shearCorrectionFactorVertical))
print()
print('Horizontal bending stiffness matrix:')
SimplifyPrint(horizontalStiffnessMatrix, factor = EIHorizontal / l**2 / (1 + shearCorrectionFactorHorizontal))
print()
print('Axial mass matrix:')
SimplifyPrint(axialMassMatrix, factor = l * linearDensity)
print()
print('Vertical bending mass matrix:')
print('Heave contribution:')
SimplifyPrint(verticalBendingHeaveContributionMassMatrix, factor = l * linearDensity / 840 / (shearCorrectionFactorVertical**2 + 2 * shearCorrectionFactorVertical + 1))
print('Pitch contribution:')
SimplifyPrint(verticalBendingPitchContributionMassMatrix, factor = (verticalMassInertiaMoment + linearDensity * (zCenterOfMass - zNeutralAxis)**2) / 30 / l / (shearCorrectionFactorVertical**2 + 2 * shearCorrectionFactorVertical + 1))
print()
print('Axial-vertical bending coupling mass matrix:')
SimplifyPrint(axialVerticalBendingCouplingMassMatrix, factor = linearDensity * (zCenterOfMass - zNeutralAxis) / (12 * (1 + shearCorrectionFactorVertical)))
print()
print('Horizontal bending mass matrix:')
print('Sway contribution')
SimplifyPrint(horizontalBendingSwayContributionMassMatrix, factor = l * linearDensity / 840 / (shearCorrectionFactorHorizontal**2 + 2 * shearCorrectionFactorHorizontal + 1))
print('Yaw contribution')
SimplifyPrint(horizontalBendingYawContributionMassMatrix, factor = horizontalMassInertiaMoment / 30 /l / (shearCorrectionFactorHorizontal**2 + 2 * shearCorrectionFactorHorizontal + 1))
print()
print('Torsion mass matrix:')
print('Bounded basis:')
SimplifyPrint(massMatrixBoundedBasis, factor = rollMassInertiaMoment + linearDensity * (zCenterOfMass - zCenterOfTwist)**2)
print('Stable basis:')
print('Free warping submatrix:')
freeWarpingSubmatrix: sym.Matrix = sym.zeros(2, 2)
freeWarpingSubmatrix[0, 0] = massMatrixStableBasis[0, 0]
freeWarpingSubmatrix[1, 0] = massMatrixStableBasis[2, 0]
freeWarpingSubmatrix[0, 1] = massMatrixStableBasis[0, 2]
freeWarpingSubmatrix[1, 1] = massMatrixStableBasis[2, 2]
SimplifyPrint(freeWarpingSubmatrix, factor = rollMassInertiaMoment + linearDensity * (zCenterOfMass - zCenterOfTwist)**2)
print('w0 w0:')
SimplifyPrint(massMatrixStableBasis[1, 1], factor = (rollMassInertiaMoment + linearDensity * (zCenterOfMass - zCenterOfTwist)**2) * l, series = True, symbol = a)
print('w0 w1:')
SimplifyPrint(massMatrixStableBasis[1, 3], factor = (rollMassInertiaMoment + linearDensity * (zCenterOfMass - zCenterOfTwist)**2) * l, series = True, symbol = a)
print('w1 w1:')
SimplifyPrint(massMatrixStableBasis[3, 3], factor = (rollMassInertiaMoment + linearDensity * (zCenterOfMass - zCenterOfTwist)**2) * l, series = True, symbol = a)
print('r0 w0:')
SimplifyPrint(massMatrixStableBasis[0, 1], factor = (rollMassInertiaMoment + linearDensity * (zCenterOfMass - zCenterOfTwist)**2) * l, series = True, symbol = a)
print('r0 w1:')
SimplifyPrint(massMatrixStableBasis[0, 3], factor = (rollMassInertiaMoment + linearDensity * (zCenterOfMass - zCenterOfTwist)**2) * l, series = True, symbol = a)
print('w0 r1:')
SimplifyPrint(massMatrixStableBasis[1, 2], factor = (rollMassInertiaMoment + linearDensity * (zCenterOfMass - zCenterOfTwist)**2) * l, series = True, symbol = a)
print('r1 w1:')
SimplifyPrint(massMatrixStableBasis[2, 3], factor = (rollMassInertiaMoment + linearDensity * (zCenterOfMass - zCenterOfTwist)**2) * l, series = True, symbol = a)
print()
print('Horizontal bending-torsion coupling mass matrix:')
print('Bounded basis:')
SimplifyPrint(horizontalBendingDispBasisTorsionBoundedBasisCouplingMassMatrix, factor = linearDensity * (zCenterOfMass - zCenterOfTwist))
print('Stable basis:')
freeWarpingCouplingSubmatrix = sym.zeros(4, 2)
freeWarpingCouplingSubmatrix[:, 0] = horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix[:, 0]
freeWarpingCouplingSubmatrix[:, 1] = horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix[:, 2]
print('Free warping submatrix:')
SimplifyPrint(freeWarpingCouplingSubmatrix, factor = linearDensity * (zCenterOfMass - zCenterOfTwist) / (1 + shearCorrectionFactorHorizontal) / 120)
print('y0 w0')
SimplifyPrint(horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix[0, 1], factor = linearDensity * l * (zCenterOfMass - zCenterOfTwist) / (1 + shearCorrectionFactorHorizontal), series = True, symbol = a)
print('psi0 w0')
SimplifyPrint(horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix[1, 1], factor = linearDensity * l**2 * (zCenterOfMass - zCenterOfTwist) / (1 + shearCorrectionFactorHorizontal), series = True, symbol = a)
print('y1 w0')
SimplifyPrint(horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix[2, 1], factor = linearDensity * l * (zCenterOfMass - zCenterOfTwist) / (1 + shearCorrectionFactorHorizontal), series = True, symbol = a)
print('psi1 w0')
SimplifyPrint(horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix[3, 1], factor = linearDensity * l**2 * (zCenterOfMass - zCenterOfTwist) / (1 + shearCorrectionFactorHorizontal), series = True, symbol = a)
print('y0 w1')
SimplifyPrint(horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix[0, 3], factor = linearDensity * l * (zCenterOfMass - zCenterOfTwist) / (1 + shearCorrectionFactorHorizontal), series = True, symbol = a)
print('psi0 w1')
SimplifyPrint(horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix[1, 3], factor = linearDensity * l**2 * (zCenterOfMass - zCenterOfTwist) / (1 + shearCorrectionFactorHorizontal), series = True, symbol = a)
print('y1 w1')
SimplifyPrint(horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix[2, 3], factor = linearDensity * l * (zCenterOfMass - zCenterOfTwist) / (1 + shearCorrectionFactorHorizontal), series = True, symbol = a)
print('psi1 w1')
SimplifyPrint(horizontalBendingDispBasisTorsionStableBasisCouplingMassMatrix[3, 3], factor = linearDensity * l**2 * (zCenterOfMass - zCenterOfTwist) / (1 + shearCorrectionFactorHorizontal), series = True, symbol = a)


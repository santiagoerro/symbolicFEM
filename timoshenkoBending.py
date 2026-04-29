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

x = sym.Symbol('x')

printType = 'python'



def PrettyPrint(expression, factor = None):
    if factor is None:
        if printType == 'latex':
            print(sym.latex(sym.simplify(expression)))
        elif printType == 'python':
            print(sym.python(sym.simplify(expression)))
        else:
            sym.pprint(sym.simplify(expression))
    else:
        if printType == 'latex':
            print(sym.latex(sym.simplify(factor)))
            print('*')
            print(sym.latex(sym.simplify(expression / factor)))
        elif printType == 'python':
            print(sym.python(sym.simplify(factor)))
            print('*')
            print(sym.python(sym.simplify(expression / factor)))
        else:
            sym.pprint(sym.simplify(factor))
            print('*')
            sym.pprint(sym.simplify(expression / factor))



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
horizontalBendingMassMatrix: sym.Matrix = linearDensity * sym.integrate(horizontalDeflectionDisplacementBasis.transpose() * horizontalDeflectionDisplacementBasis, (x, 0, l))


twistingBoundedBasis: sym.Matrix = sym.Matrix([1, x/l, sym.exp(a * (x/l - 1)), sym.exp(-a * x/l)]).transpose()

boundedBasisCoefsToConditionsMatrix: sym.Matrix = sym.zeros(4, 4)
boundedBasisCoefsToConditionsMatrix[0, :] = twistingBoundedBasis.subs(x, 0)
boundedBasisCoefsToConditionsMatrix[1, :] = twistingBoundedBasis.subs(x, l)
boundedBasisCoefsToConditionsMatrix[2, :] = sym.diff(twistingBoundedBasis, x).subs(x, 0) * l
boundedBasisCoefsToConditionsMatrix[3, :] = sym.diff(twistingBoundedBasis, x).subs(x, l) * l
stableBasisToBoundedBasisCoefsMatrix: sym.Matrix = sym.zeros(4,4)
stableBasisToBoundedBasisCoefsMatrix[:, 0] = sym.Matrix([1, -1, 0, 0])
stableBasisToBoundedBasisCoefsMatrix[:, 1] = sym.Matrix([0, 1, 0, 0])
stableBasisToBoundedBasisCoefsMatrix[:, 2] = sym.simplify(boundedBasisCoefsToConditionsMatrix.LDLsolve(sym.Matrix([0, 0, 1, 0])))
stableBasisToBoundedBasisCoefsMatrix[:, 3] = sym.simplify(boundedBasisCoefsToConditionsMatrix.LDLsolve(sym.Matrix([0, 0, 0, 1])))

# twistingStableBasis: sym.Matrix = twistingBoundedBasis * stableBasisToBoundedBasisCoefsMatrix
# PrettyPrint(twistingStableBasis)

massMatrixBoundedBasis = sym.integrate(twistingBoundedBasis.transpose() * twistingBoundedBasis, (x, 0, l))
massMatrixStableBasis = sym.simplify(stableBasisToBoundedBasisCoefsMatrix.transpose() * massMatrixBoundedBasis * stableBasisToBoundedBasisCoefsMatrix)


bendingverticalDeflectionDisplacementBasisTwistingBoundedBasis = sym.Matrix.hstack(verticalDeflectionDisplacementBasis, twistingBoundedBasis)

bendingTorsionCoupledMassMatrixDisplacementBoundedBasis = sym.integrate(bendingverticalDeflectionDisplacementBasisTwistingBoundedBasis.transpose() * bendingverticalDeflectionDisplacementBasisTwistingBoundedBasis, (x, 0, l))



print()
print('Axial force:')
PrettyPrint(axialForceDisplacementBasis, factor = EA/l)
print()
print('Axial stiffness matrix:')
PrettyPrint(axialStiffnessMatrix, factor = EA/l)
print()
print('Vertical bending moment:')
PrettyPrint(verticalBendingMomentDisplacementBasis, factor = EIVertical / l**2 / (1 + shearCorrectionFactorVertical))
print()
print('Horizontal bending moment:')
PrettyPrint(horizontalBendingMomentDisplacementBasis, factor = EIHorizontal / l**2 / (1 + shearCorrectionFactorHorizontal))
print()
print('Vertical shear force:')
PrettyPrint(verticalShearForceDisplacementBasis, factor = EIVertical / l**3 / (1 + shearCorrectionFactorVertical))
print()
print('Horizontal shear force:')
PrettyPrint(horizontalShearForceDisplacementBasis, factor = EIHorizontal / l**3 / (1 + shearCorrectionFactorHorizontal))
print()
print('Vertical bending stiffness matrix:')
PrettyPrint(verticalStiffnessMatrix, factor = EIVertical / l**2 / (1 + shearCorrectionFactorVertical))
print()
print('Horizontal bending stiffness matrix:')
PrettyPrint(horizontalStiffnessMatrix, factor = EIHorizontal / l**2 / (1 + shearCorrectionFactorHorizontal))
print()
print('Axial mass matrix:')
PrettyPrint(axialMassMatrix, factor = l * linearDensity)
print()
print('Vertical bending mass matrix:')
print('Heave contribution:')
PrettyPrint(verticalBendingHeaveContributionMassMatrix, factor = l * linearDensity / 840 / (shearCorrectionFactorVertical**2 + 2 * shearCorrectionFactorVertical + 1))
print('Pitch contribution:')
PrettyPrint(verticalBendingPitchContributionMassMatrix, factor = (verticalMassInertiaMoment + linearDensity * (zCenterOfMass - zNeutralAxis)**2) / 30 / l / (shearCorrectionFactorVertical**2 + 2 * shearCorrectionFactorVertical + 1))
print()
print('Axial-vertical bending coupling mass matrix:')
PrettyPrint(axialVerticalBendingCouplingMassMatrix, factor = linearDensity * (zCenterOfMass - zNeutralAxis) / (12 * (1 + shearCorrectionFactorVertical)))
print()
print('Horizontal bending mass matrix:')
PrettyPrint(horizontalBendingMassMatrix, factor = l * linearDensity / 840 / (shearCorrectionFactorHorizontal**2 + 2 * shearCorrectionFactorHorizontal + 1))
print()
print('Torsion mass matrix:')
print('Bounded basis:')
PrettyPrint(massMatrixBoundedBasis)
print('Stable basis:')
print('Free warping submatrix:')
PrettyPrint(massMatrixStableBasis[0:2, 0:2])
print('w0 w0:')
PrettyPrint(massMatrixStableBasis[2, 2])
PrettyPrint(sym.series(massMatrixStableBasis[2, 2], a, 0, 10))
print('w0 w1:')
PrettyPrint(massMatrixStableBasis[2, 3])
PrettyPrint(sym.series(massMatrixStableBasis[2, 3], a, 0, 10))
print('w1 w1:')
PrettyPrint(massMatrixStableBasis[3, 3])
PrettyPrint(sym.series(massMatrixStableBasis[3, 3], a, 0, 10))
print('r0 w0:')
PrettyPrint(massMatrixStableBasis[0, 2])
PrettyPrint(sym.series(massMatrixStableBasis[0, 2], a, 0, 10))
print('r0 w1:')
PrettyPrint(massMatrixStableBasis[0, 3])
PrettyPrint(sym.series(massMatrixStableBasis[0, 3], a, 0, 10))
print('r1 w0:')
PrettyPrint(massMatrixStableBasis[1, 2])
PrettyPrint(sym.series(massMatrixStableBasis[1, 2], a, 0, 10))
print('r1 w1:')
PrettyPrint(massMatrixStableBasis[1, 3])
PrettyPrint(sym.series(massMatrixStableBasis[1, 3], a, 0, 10))

# print()
# print('Bending-torsion coupled mass matrix:')
# PrettyPrint(bendingTorsionCoupledMassMatrixDisplacementBoundedBasis.subs(delta, 0))


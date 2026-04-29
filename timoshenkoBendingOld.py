import sympy as sym


l = sym.Symbol('l', positive = True)
delta = sym.Symbol('delta', positive = True)
a = sym.Symbol('a', positive = True)
EI = sym.Symbol('EI')
linearDensity = sym.Symbol('lambda')

x = sym.Symbol('x')


deflectionMonomialBasis: sym.Matrix = sym.Matrix([1, x, x**2, x**3]).transpose()
rotationMonomialBasis: sym.Matrix = sym.diff(deflectionMonomialBasis, x) + delta * l**2 / 12 * sym.diff(deflectionMonomialBasis, x, 3)

monomialBasisCoefsToDisplacementsMatrix: sym.Matrix = sym.zeros(4, 4)
monomialBasisCoefsToDisplacementsMatrix[0, :] = deflectionMonomialBasis.subs(x, 0)
monomialBasisCoefsToDisplacementsMatrix[1, :] = rotationMonomialBasis.subs(x, 0)
monomialBasisCoefsToDisplacementsMatrix[2, :] = deflectionMonomialBasis.subs(x, l)
monomialBasisCoefsToDisplacementsMatrix[3, :] = rotationMonomialBasis.subs(x, l)

displacementsToMonomialBasisCoefsMatrix: sym.Matrix = monomialBasisCoefsToDisplacementsMatrix.inv()
deflectionDisplacementBasis: sym.Matrix = deflectionMonomialBasis * displacementsToMonomialBasisCoefsMatrix
rotationDisplacementBasis: sym.Matrix = rotationMonomialBasis * displacementsToMonomialBasisCoefsMatrix

bendingMomentDisplacementBasis: sym.Matrix = EI * sym.diff(rotationDisplacementBasis, x)
shearForceDisplacementBasis: sym.Matrix = 12 * EI / (delta * l**2) * (sym.diff(deflectionDisplacementBasis, x) - rotationDisplacementBasis)

stiffnessMatrix: sym.Matrix = sym.zeros(4, 4)
stiffnessMatrix[0, :] = -shearForceDisplacementBasis.subs(x, 0)
stiffnessMatrix[1, :] = -bendingMomentDisplacementBasis.subs(x, 0)
stiffnessMatrix[2, :] = shearForceDisplacementBasis.subs(x, l)
stiffnessMatrix[3, :] = bendingMomentDisplacementBasis.subs(x, l)

massMatrix: sym.Matrix = linearDensity * sym.integrate(deflectionDisplacementBasis.transpose() * deflectionDisplacementBasis, (x, 0, l))


twistingBoundedBasis: sym.Matrix = sym.Matrix([1 - x/l, x/l, sym.exp(a * (x/l - 1)), sym.exp(-a * x/l)]).transpose()

massMatrixBoundedBasis = sym.integrate(twistingBoundedBasis.transpose() * twistingBoundedBasis, (x, 0, l))


bendingDeflectionDisplacementBasisTwistingBoundedBasis = sym.Matrix.hstack(deflectionDisplacementBasis, twistingBoundedBasis)

bendingTorsionCoupledMassMatrixDisplacementBoundedBasis = sym.integrate(bendingDeflectionDisplacementBasisTwistingBoundedBasis.transpose() * bendingDeflectionDisplacementBasisTwistingBoundedBasis, (x, 0, l))



def PrettyPrint(expression, factor = None):
    if factor is None:
        sym.pprint(sym.simplify(expression))
    else:
        sym.pprint(sym.simplify(factor))
        print('*')
        sym.pprint(sym.simplify(expression / factor))



print()
print('Bending moment:')
PrettyPrint(bendingMomentDisplacementBasis, factor = EI / l**2 / (1 + delta))
print()
print('Shear force:')
PrettyPrint(shearForceDisplacementBasis, factor = EI / l**3 / (1 + delta))
print()
print('Bending stiffness matrix:')
PrettyPrint(stiffnessMatrix, factor = EI / l**2 / (1 + delta))
print()
print('Bending mass matrix:')
PrettyPrint(massMatrix, factor = l * linearDensity / 840 / (delta**2 + 2 * delta + 1))
print()
print()
print('Torsion mass matrix:')
PrettyPrint(massMatrixBoundedBasis)
print()
print()
print('Bending-torsion coupled mass matrix:')
PrettyPrint(bendingTorsionCoupledMassMatrixDisplacementBoundedBasis.subs(delta, 0))


import sympy as sp
import sympy.physics.units.quantities as sq
from sympy.physics.quantum.constants import hbar
import copy

PSI_FUNCTION = sp.Function("psi")
POTENTIAL_FUNCTION = sp.Function("V")
TOTAL_ENERGY_SYMBOL = sq.Symbol('E', nonzero = True, real = True, positive = True)
MASS_SYMBOL = sq.Quantity('m', positive = True, nonzero = True)
POSITION_SYMBOL = sp.Symbol('x', positive = True)
DEFAULT_CONJUGATE_NOT_SQUARED_ABSOLUTE_VALUE = True
DEFAULT_NORMALIZATION_VALUE = 1
INDEFINITE_NORMALIZATION_INTEGRAL = None
DEFAULT_CONJUGATE_NOT_SQUARED_ABSOLUTE_VALUE = True

def makeDistance(region): 
    return sp.Symbol(f'L_{region}', real = True, finite = True)

def makeStartDistance(region): 
    return 0 if region == 0 else sp.Symbol(f'L_{region - 1}', positive = True, real = True)

class RegionSymbols: 
    def __init__(self, regionNumber): 
        self.regionNumber = regionNumber
        self.waveFunction = sp.Function(f'\psi_{self.regionNumber}')
        self.constants = [
                sp.Symbol('C_{{' + str(self.regionNumber) + '}_t}'), 
                sp.Symbol('C_{{' + str(self.regionNumber) + '}_r}')
            ]
        self.distance = makeDistance(self.regionNumber)
        self.harmonicConstant = sp.Symbol(f'k_{self.regionNumber}', real = True)
        self.boundry = sp.Symbol(f'B_{self.regionNumber}')
        self.normalizationConstant = sp.Symbol(f'N_{self.regionNumber}')
        self.partSummeries = [
                sp.Function('\psi_{{' + str(self.regionNumber) + '}_t}'), 
                sp.Function('\psi_{{' + str(self.regionNumber) + '}_r}')
            ]
        self.mass = sq.Quantity('m', positive = True, nonzero = True)
        self.constantPotential = sp.Symbol(f'V_{self.regionNumber}', positive = True, real = True)
        self.startDistance = makeStartDistance(self.regionNumber)

def makeCoefficents(first, second): 
    transmission = (sp.Abs(second.constants[0]) ** 2) / (sp.Abs(first.constants[0]) ** 2)
    reflection = (sp.Abs(first.constants[1]) ** 2) / (sp.Abs(first.constants[0]) ** 2)
    return transmission, reflection

def makeCoefficentsFromHarmonicConstants(first, second): 
    transmission = 2 * first.harmonicConstant / (first.harmonicConstant + second.harmonicConstant)
    reflection =  (first.harmonicConstant - second.harmonicConstant) / (first.harmonicConstant + second.harmonicConstant)
    return transmission, reflection

def removeAbsAndConjugate(toConjugate): 
    wild = sp.Wild('wild')
    return toConjugate.replace(sp.Abs(wild), sp.conjugate(wild))

def removeAbs(toRemoveFrom): 
    wild = sp.Wild('wild')
    return toRemoveFrom.replace(sp.Abs(wild), wild)

def ratioFromNormalizeRatio(ratio): 
    ratio_ = sp.sqrt(ratio)
    return removeAbs(ratio), removeAbsAndConjugate(ratio)

def makeScatteringMatrix(transmission, reflection): 
    transmission, transmissionConjugate = ratioFromNormalizeRatio(sp.sqrt(transmission))
    reflection, reflectionConjugate = ratioFromNormalizeRatio(sp.sqrt(reflection))
    return sp.Matrix([
            [reflectionConjugate, transmission], 
            [transmissionConjugate, reflection]
        ])

def makeTransmissionMatrix(scatteringMatrix): 
    return sp.Matrix([
            [
                    scatteringMatrix.col(1).row(0) - (scatteringMatrix.col(1).row(1) 
                            * (scatteringMatrix.col(0).row(0)) * (scatteringMatrix.col(0).row(1) ** (-1))), 
                    scatteringMatrix.col(1).row(1) * (scatteringMatrix.col(0).row(1) ** (-1))
            ], 
            [
                    -scatteringMatrix.col(0).row(0) * (scatteringMatrix.col(0).row(1) ** (-1)), 
                    scatteringMatrix.col(0).row(1) ** (-1)
            ]
        ])

def makeBarrierMatrix(from_ : RegionSymbols): 
    return sp.Matrix([
            [from_.partSummeries[0](from_.startDistance)], 
            [from_.partSummeries[1](from_.startDistance)]
        ])

def transfer(from_, to): 
    transmission, reflection = makeCoefficentsFromHarmonicConstants(from_, to)
    scatteringMatrix = makeScatteringMatrix(transmission, reflection)
    transferMatrix = makeTransmissionMatrix(scatteringMatrix)
    return sp.Eq(makeBarrierMatrix(from_), transferMatrix * makeBarrierMatrix(to))

def parameterizeTransferForSimulation(symbols, transfer, baseName): 
    inputs = [
            sp.Symbol(baseName + "_t"), 
            sp.Symbol(baseName + "_r")
        ]
    transfer = transfer.replace(symbols.partSummeries[0](symbols.startDistance), inputs[0])
    transfer = transfer.replace(symbols.partSummeries[1](symbols.startDistance), inputs[1])
    return inputs, transfer

def lambdifyTransfer(transferSymbols, harmonicConstants): 
    display(transferSymbols[0][0])
    display(transferSymbols[1].rhs.row(0)[0])
    display(transferSymbols[0][1])
    display(transferSymbols[1].rhs.row(1)[0])
    return [
        sp.lambdify([*tuple(transferSymbols[0]), *tuple(harmonicConstants)], transferSymbols[1].rhs.row(0)[0]), 
        sp.lambdify([*tuple(transferSymbols[0]), *tuple(harmonicConstants)], transferSymbols[1].rhs.row(1)[0])
    ]

def calculateTransfers(regionFunctions, initialValues, log, harmonicConstants): 
    if len(regionFunctions) == 1:
        log.append((
                regionFunctions[0][0](*tuple(initialValues), *tuple(harmonicConstants)), \
                regionFunctions[0][1](*tuple(initialValues), *tuple(harmonicConstants))
            ))
        return log[-1]
    elif len(regionFunctions) > 0: 
        transmission, reflection = calculateTransfers(regionFunctions[1:], initialValues, log, harmonicConstants)
        log.append((
                regionFunctions[0][0](transmission, reflection, *tuple(harmonicConstants)), \
                regionFunctions[0][1](transmission, reflection, *tuple(harmonicConstants))
            ))
        return log[-1]
    else: 
        assert True, "Please input at least one pair of region functions"
        
def constantPotentialTimeIndependentSchroedingerEquation1D( 
            regionSymbols : RegionSymbols, 
            totalEnergy = TOTAL_ENERGY_SYMBOL, 
            reducedPlanckConstant = hbar, 
            position = POSITION_SYMBOL 
        ): 
    return sp.Eq( 
            ( ( -( reducedPlanckConstant ** 2 ) / ( 2 * regionSymbols.mass ) ) \
                    * sp.Derivative(regionSymbols.waveFunction(position), (position, 2)))
                    + (regionSymbols.constantPotential * regionSymbols.waveFunction(position)), 
            totalEnergy * regionSymbols.waveFunction(position)
        )
        
def simpleWaveFunctionNormalization( 
            from_, toOrIndefinite, 
            regionSymbols : RegionSymbols, 
            position = POSITION_SYMBOL, 
            conjugateNotSquaredAbsoluteValue 
                    = DEFAULT_CONJUGATE_NOT_SQUARED_ABSOLUTE_VALUE 
        ): 
    integralFunction = regionSymbols.waveFunction(position) * sp.conjugate(regionSymbols.waveFunction(position)) \
            if conjugateNotSquaredAbsoluteValue \
            else sp.Abs(regionSymbols.waveFunction(position)) ** 2
    if toOrIndefinite: 
        return sp.Eq(
                sp.Integral(
                        integralFunction, 
                        (position, from_, toOrIndefinite)
                    ),
                regionSymbols.normalizationConstant
            )
    else:
        return sp.Eq(
                sp.Integral(integralFunction), 
                regionSymbols.normalizationConstant
            )

def manipulate(equation, operation): 
    return sp.Eq(operation(equation.lhs), operation(equation.rhs))

def secondDerivative(regionSymbols : RegionSymbols, position = POSITION_SYMBOL): 
    return sp.Derivative(regionSymbols.waveFunction(position), (position, 2))

def constantFactor(
            regionSymbols : RegionSymbols, 
            mass = MASS_SYMBOL, 
            reducedPlanckConstant = hbar
        ):
    return (reducedPlanckConstant ** 2) / (2 * mass)

def secondDerivativeTerm(
            regionSymbols : RegionSymbols, 
            position = POSITION_SYMBOL, 
            mass = MASS_SYMBOL, 
            reducedPlanckConstant = hbar
        ):
    return constantFactor(regionSymbols) * secondDerivative(regionSymbols)

def extractHarmonicConstant(regionSymbols, equation): 
    waveFunction = regionSymbols.waveFunction(POSITION_SYMBOL)
    solveForWaveFunction = sp.Eq(waveFunction, sp.solve(equation, regionSymbols.waveFunction(POSITION_SYMBOL))[0])
    return sp.Eq(
            regionSymbols.harmonicConstant, 
            manipulate(solveForWaveFunction, lambda step : sp.sqrt((step / secondDerivative(regionSymbols)) ** (-1))).refine().rhs
        )


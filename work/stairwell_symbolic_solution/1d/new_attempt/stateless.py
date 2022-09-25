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

def makeFromToSymbol(baseName, index): 
    return {
            "next" : sp.Symbol(f"{baseName}_{index}_{index + 1}"), 
            "previous" : sp.Symbol(f"{baseName}_{index}_{index - 1}")
        }

class NextPrevious: 
    def __init__(self, transmissionCoefficent, reflectionCoefficient): 
        self.transmissionCoefficent = transmissionCoefficent
        self.reflectionCoefficient = reflectionCoefficient

class RegionSymbols: 
    def __init__(self, regionIndex): 
        self.regionIndex = regionIndex
        self.waveFunction = sp.Function(f'\psi_{self.regionIndex}')
        self.constants = [
                sp.Symbol('C_{{' + str(self.regionIndex) + '}_t}'), 
                sp.Symbol('C_{{' + str(self.regionIndex) + '}_r}')
            ]
        self.distance = makeDistance(self.regionIndex)
        self.harmonicConstant = sp.Symbol(f'k_{self.regionIndex}', real = True)
        self.boundry = sp.Symbol(f'B_{self.regionIndex}')
        self.normalizationConstant = sp.Symbol(f'N_{self.regionIndex}')
        self.partSummeries = [
                sp.Function('\psi_{{' + str(self.regionIndex) + '}_t}'), 
                sp.Function('\psi_{{' + str(self.regionIndex) + '}_r}')
            ]
        self.mass = sq.Quantity('m', positive = True, nonzero = True)
        self.constantPotential = sp.Symbol(f'V_{self.regionIndex}', positive = True, real = True)
        self.startDistance = makeStartDistance(self.regionIndex)
        self.transmissionCoefficents = makeFromToSymbol("T", self.regionIndex)
        self.reflectionCoefficents = makeFromToSymbol("R", self.regionIndex)
        self.next = NextPrevious(self.transmissionCoefficents["next"], self.reflectionCoefficents["next"])
        self.previous = NextPrevious(self.transmissionCoefficents["previous"], self.reflectionCoefficents["previous"])

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
    noAbs = removeAbs(ratio)
    return ratio, sp.conjugate(ratio)

def makeScatteringMatrix(from_, to): 
    ratios = from_.previous if from_.regionIndex > to.regionIndex else from_.next
    transmission, transmissionConjugate = ratioFromNormalizeRatio(ratios.transmissionCoefficent)
    reflection, reflectionConjugate = ratioFromNormalizeRatio(ratios.reflectionCoefficient)
    return sp.Matrix([
            [sp.sqrt(reflectionConjugate), sp.sqrt(transmission)], 
            [sp.sqrt(transmissionConjugate), sp.sqrt(reflection)]
        ])

def makeTransferMatrix(scatteringMatrix): 
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

def transferWithHarmonicConstants(from_, to): 
    transmission, reflection = makeCoefficentsFromHarmonicConstants(from_, to)
    scatteringMatrix = makeScatteringMatrix(transmission, reflection)
    return generalTransferFromScatteringMatrix(from_, to, scatteringMatrix)

def generalTransfer(from_, to): 
    return generalTransferFromScatteringMatrix(from_, to, makeScatteringMatrix(from_, to))

def generalTransferFromScatteringMatrix(from_, to, scatteringMatrix): 
    transferMatrix = makeTransferMatrix(scatteringMatrix)
    return {
            "from" : from_, 
            "to" : to, 
            "scatteringMatrix" : scatteringMatrix, 
            "transferMatrix" : transferMatrix, 
            "transfer" : sp.Eq(makeBarrierMatrix(from_), transferMatrix * makeBarrierMatrix(to))
        }

def parameterizeTransferForSimulation(transferData): 
    inputs = [
            sp.Symbol(str(symbols.distance) + "_t"), 
            sp.Symbol(str(symbols.distance) + "_r")
        ]
    transfer = transfer.replace(symbols.partSummeries[0](symbols.startDistance), inputs[0])
    transfer = transfer.replace(symbols.partSummeries[1](symbols.startDistance), inputs[1])
    return { "inputs" : inputs, "paramterizedTransfer" : transfer }

def lambdifyTransfer(transferSymbols, scatteringParameters): 
    parameters = [*tuple(transferSymbols["inputs"]), *tuple(scatteringParameters)]
    calculation = transferSymbols["transfer"].rhs
    result = transferSymbols["transfer"].lhs
    return {
            "parameters" : parameters, 
            str(result.row(0)[0]) + "Function" : sp.lambdify(parameters, calculation.row(0)[0]), 
            str(result.row(1)[0]) + "Function" : sp.lambdify(parameters, calculation.row(1)[0])
        }

#def calculateTransfers(transferFunctions, initialValues, scatteringConstants): 
    

def calculateTransfers(regionFunctions, initialValues, log, harmonicConstants, regionBoundrySymbols): 
    if len(regionFunctions) == 1:
        log.append((
                regionFunctions[0][0](*tuple(initialValues), *tuple(harmonicConstants)), \
                regionFunctions[0][1](*tuple(initialValues), *tuple(harmonicConstants)), \
                regionBoundrySymbols[0][0]
            ))
        return log[-1]
    elif len(regionFunctions) > 0: 
        transmission, reflection, nextRegionBoundrySymbols = calculateTransfers(
                regionFunctions[1:], 
                initialValues, 
                log, 
                harmonicConstants, 
                regionBoundrySymbols[1:]
            )
        log.append((
                regionFunctions[0][0](transmission, reflection, *tuple(harmonicConstants)), \
                regionFunctions[0][1](transmission, reflection, *tuple(harmonicConstants)), \
                regionBoundrySymbols[0][0]
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

#def solveODE(symbols : RegionSymbol, schrodingerEquation):
#        


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


import sympy as sp
import sympy.physics.units.quantities as sq
from sympy.physics.quantum.constants import hbar
import copy
from collections import OrderedDict
import inspect

PSI_FUNCTION = sp.Function("psi")
POTENTIAL_FUNCTION = sp.Function("V")
TOTAL_ENERGY_SYMBOL = sq.Symbol('E_{total}', nonzero = True, real = True, positive = True)
MASS_SYMBOL = sq.Symbol('m', positive = True, nonzero = True)
POSITION_SYMBOL = sp.Symbol('x', positive = True)
DEFAULT_CONJUGATE_NOT_SQUARED_ABSOLUTE_VALUE = True
DEFAULT_NORMALIZATION_VALUE = 1
INDEFINITE_NORMALIZATION_INTEGRAL = None
DEFAULT_CONJUGATE_NOT_SQUARED_ABSOLUTE_VALUE = True

RADICAND_TAG = "radicand"
INVERSE_COEFFICIENT_TAG = "inverseCoefficent"
OFFSET_TAG = "offset"
SCALAR_TAG = "scalar"

def makeDistance(region): 
    return sp.Symbol(f'L_{region}', real = True, finite = True)

def makeStartDistance(region): 
    return sp.Symbol(f'L_{region - 1}', positive = True, real = True)

def makeFromToSymbol(baseName, index): 
    return {
            "next" : sp.Symbol(f"{baseName}_{index}_{index + 1}"), 
            "previous" : sp.Symbol(f"{baseName}_{index}_{index - 1}")
        }

class RegionCoefficients: 
    def __init__(self, transmissionCoefficent, reflectionCoefficient): 
        self.transmissionCoefficent = transmissionCoefficent
        self.reflectionCoefficient = reflectionCoefficient

def isIdentifier(identifier):
    firstCharacter = ord(identifier[0])
    lowercase = (firstCharacter >= ord('a') and firstCharacter <= ord('z'))
    uppercase = (firstCharacter >= ord('A') and firstCharacter <= ord('Z'))
    traditionalIdentifier = lowercase == True or uppercase == True or firstCharacter == ord('_')
    if identifier.startswith("conjugate") or identifier.startswith("exp("):
        return False
    return traditionalIdentifier

def functionNameFromIdentifier(identifier): 
    return identifier + "Function"

def orderNames(names, sortIndex = 0): 
    sameStart = {}
    results = []
    for name in names: 
        if sortIndex < len(name): 
            sortCharacter = name[sortIndex]
            if sortCharacter in sameStart.keys(): 
                sameStart[sortCharacter].append(name)
            else: 
                sameStart[sortCharacter] = [name]
        else: 
            results.append(name)
    for ordered in [orderNames(sameStart[character], sortIndex + 1) for character in sorted(list(sameStart))]: 
        results += ordered
    return results 

def symbolicToIdentifier(symbolic):
    identifier = str(symbolic)
    replacementList = {
            '{' : "OCB", 
            '}' : "CCB", 
            '(' : "OP", 
            ')' : "CP", 
            '[' : "OSB", 
            ']' : "CSB", 
            '-' : "N", 
            "Backslash" : "\\"
        }
    for toReplace, replacement in replacementList.items(): 
        identifier = identifier.replace(toReplace, replacement)
    return identifier

def substituteIdentifierAtomsList(symbolic):
    substitutionList = {}
    atoms = symbolic.atoms()
    atoms |= symbolic.atoms(sp.Function)
    for atom in atoms: 
        if isIdentifier(str(atom)) == True: 
            substitutionList[atom] = symbolicToIdentifier(atom)
    return substitutionList

def substituteIdentifierAtoms(symbolic):
    return symbolic.subs(substituteIdentifierAtomsList(symbolic))

def makeBoundrySymbol(regionIndex : int): 
    return sp.Symbol(f'B_{regionIndex}')

class RegionSymbols: 
    def __init__(self, regionIndex): 
        self.regionIndex = regionIndex
        self.waveFunction = sp.Function(f'psi_{self.regionIndex}')
        self.constants = [
                sp.Symbol('C_{{' + str(self.regionIndex) + '}_t}'), 
                sp.Symbol('C_{{' + str(self.regionIndex) + '}_r}')
            ]
        self.distance = makeDistance(self.regionIndex)
        self.harmonicConstant = sp.Symbol(f'k_{self.regionIndex}', real = True)
        self.boundry = makeBoundrySymbol(self.regionIndex)
        self.nextBoundry = makeBoundrySymbol(self.regionIndex + 1)
        self.normalizationConstant = sp.Symbol(f'N_{self.regionIndex}')
        self.partSummeries = [
                sp.Function(f'psi_{str(self.regionIndex)}_t'), 
                sp.Function(f'psi_{str(self.regionIndex)}_r')
            ]
        self.mass = MASS_SYMBOL
        self.constantPotential = sp.Symbol(f'V_{self.regionIndex}', positive = True, real = True)
        self.startDistance = makeStartDistance(self.regionIndex)
        self.transmissionCoefficents = makeFromToSymbol("T", self.regionIndex)
        self.reflectionCoefficents = makeFromToSymbol("R", self.regionIndex)
        self.next = RegionCoefficients(self.transmissionCoefficents["next"], self.reflectionCoefficents["next"])
        self.previous = RegionCoefficients(self.transmissionCoefficents["previous"], self.reflectionCoefficents["previous"])

def makeCoefficents(first, second): 
    transmission = (sp.Abs(second.constants[0]) ** 2) / (sp.Abs(first.constants[0]) ** 2)
    reflection = (sp.Abs(first.constants[1]) ** 2) / (sp.Abs(first.constants[0]) ** 2)
    return RegionCoefficients(transmission, reflection)

def makeCoefficentsFromHarmonicConstants(first, second): 
    transmission = 2 * first.harmonicConstant / (first.harmonicConstant + second.harmonicConstant)
    reflection =  (first.harmonicConstant - second.harmonicConstant) / (first.harmonicConstant + second.harmonicConstant)
    return RegionCoefficients(transmission, reflection)

def removeAbsAndConjugate(toConjugate): 
    wild = sp.Wild('wild')
    return toConjugate.replace(sp.Abs(wild), sp.conjugate(wild))

def removeAbs(toRemoveFrom): 
    wild = sp.Wild('wild')
    return toRemoveFrom.replace(sp.Abs(wild), wild)

def ratioFromNormalizeRatio(ratio): 
    noAbs = removeAbs(ratio)
    return ratio, sp.conjugate(ratio)

def selectRegionCoefficients(
            from_ : RegionSymbols, 
            to : RegionSymbols, 
            coefficents : RegionCoefficients = None
        ) -> RegionCoefficients: 
    return coefficents if coefficents else (from_.previous if from_.regionIndex > to.regionIndex else from_.next)

def makeScatteringMatrix(from_, to, coefficents : RegionCoefficients = None): 
    ratios = selectRegionCoefficients(from_, to, coefficents)
    transmission, transmissionConjugate = ratioFromNormalizeRatio(ratios.transmissionCoefficent)
    reflection, reflectionConjugate = ratioFromNormalizeRatio(ratios.reflectionCoefficient)
    return {
            "inputs" : [
                    symbolicToIdentifier(ratios.transmissionCoefficent), 
                    symbolicToIdentifier(ratios.reflectionCoefficient)
                ], 
            "scatteringMatrix" : sp.Matrix([
                    [sp.sqrt(reflectionConjugate), sp.sqrt(transmission)], 
                    [sp.sqrt(transmissionConjugate), sp.sqrt(reflection)]
                ])
        }

def makeTransferMatrix(scatteringMatrix_): 
    scatteringMatrix = scatteringMatrix_["scatteringMatrix"]
    scatteringMatrix_["transferMatrix"] = sp.Matrix([
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
    return scatteringMatrix_

def makeBarrierMatrix(from_ : RegionSymbols): 
    transmission = from_.partSummeries[0](from_.startDistance)
    reflection = from_.partSummeries[1](from_.startDistance)
    return {
            "inputs" : [
                    symbolicToIdentifier(transmission), 
                    symbolicToIdentifier(reflection)
                ], 
            "matrix" : sp.Matrix([
                    [transmission], 
                    [reflection]
                ])
        }

def transferWithHarmonicConstants(from_, to): 
    transmission, reflection = makeCoefficentsFromHarmonicConstants(from_, to)
    scatteringMatrix = makeScatteringMatrix(transmission, reflection)
    return generalTransferFromScatteringMatrix(from_, to, scatteringMatrix)

def generalTransfer(from_, to): 
    return generalTransferFromScatteringMatrix(from_, to, makeScatteringMatrix(from_, to))

def generalTransferFromScatteringMatrix(from_, to, scatteringMatrix): 
    transferMatrix = makeTransferMatrix(scatteringMatrix)
    input_ = makeBarrierMatrix(to)
    result = makeBarrierMatrix(from_)
    transferMatrix["inputs"] += result["inputs"]
    transferMatrix["outputs"] = input_["inputs"]
    return {
            "from" : from_, 
            "to" : to, 
            "matricies" : transferMatrix, 
            "transfer" : sp.Eq(input_["matrix"], transferMatrix["transferMatrix"] * result["matrix"])
        }

def lambdifyTransfer(transferData): 
    parameters = transferData["matricies"]["inputs"]
    outputs = transferData["matricies"]["outputs"]
    calculation = transferData["transfer"].rhs
    result = transferData["transfer"].lhs
    transferData["functionData"] = { 
            "paramters" : parameters, 
            "functions" : {
                    functionNameFromIdentifier(outputs[0]) : sp.lambdify(
                            parameters, 
                            substituteIdentifierAtoms(calculation.row(0)[0])
                        ), 
                    functionNameFromIdentifier(outputs[1]) : sp.lambdify(
                            parameters, 
                            substituteIdentifierAtoms(calculation.row(1)[0])
                        )
                }
        }
    return transferData

def inputCoefficients(functions, inputs, coefficients): 
    assert len(functions) == 2, ("Wrong number of of functions to predict the wave function (reflective and transmission)"
                                + "at a certain coordinate, are you using more than 2 dimensions?")
    satisfiedArguments, waveFunctionInputs = satisfyParameterDict(inputs, coefficients).values()
    assert len(waveFunctionInputs) == 2, "Insufficiant arguemnts for coefficents"
    return satisfiedArguments, waveFunctionInputs

def performTransfersImplementation(
            transfers : list[dict], 
            previousTransferValues, 
            transferValues, 
            transmissionReflectionCoefficients : dict
        ) -> dict: 
    functionData = transfers[0]['functionData']
    inputs = transfers[0]['matricies']['inputs']
    outputs = transfers[0]['matricies']['outputs']
    functions = functionData['functions']
    satisfiedArguments, waveFunctionInputs = inputCoefficients(
            functions, 
            inputs, 
            transmissionReflectionCoefficients[0]
        )
    arguments, remainingArguments = satisfyParameterDict(waveFunctionInputs, previousTransferValues).values()
    assert len(remainingArguments) == 0, "Parameter not satisfied!"
    arguments = tuple(satisfiedArguments.values()) + tuple(arguments.values())
    waveFunctions = tuple(functions.keys())
    transferValues |= previousTransferValues
    thisTransferValues = {
            outputs[0] : functions[waveFunctions[0]](*arguments), 
            outputs[1] : functions[waveFunctions[1]](*arguments)
        }
    if len(transfers) <= 1:
        return transferValues | thisTransferValues
    else: 
        return performTransfersImplementation(
                transfers[1:], 
                thisTransferValues, 
                transferValues, 
                transmissionReflectionCoefficients[1:]
            )


def performTransfers(
            transfers : list[dict], 
            initialTransmission, 
            initialReflection, 
            transmissionReflectionCoefficients : dict
        ) -> dict: 
    functions = transfers[0]['functionData']['functions']
    inputs = transfers[0]['matricies']['inputs']
    satisfiedArguments, waveFunctionInputs = inputCoefficients(
            functions, 
            inputs, 
            transmissionReflectionCoefficients[0]
        )
    return performTransfersImplementation(
            transfers, 
            {
                    waveFunctionInputs[0] : initialTransmission, 
                    waveFunctionInputs[1] : initialReflection
            }, 
            {}, 
            transmissionReflectionCoefficients
        )

def satisfyParameterDict(parameters, inputMapping): 
    arguments = {}
    isMappingType = type(inputMapping) is dict or type(inputMapping) is OrderedDict
    assert isMappingType, "inputMapping is not of type dict, need to map names of paramters to values"
    for parameter, argument in inputMapping.items(): 
        satisfied = parameter == parameters[0]
        parameters = parameters[1:]
        if satisfied == False: 
            assert satisfied, ("satisfyParameters: Parameter not in parameter list or is not "
                    "the current parameter (parameter order must be preserved)!")
        else: 
            arguments[parameter] = argument
    return {"arguments" : arguments, "remainingParameters" : parameters}

def defaultTransmissionReflectionCoefficientGenerator(
            from_ : RegionSymbols, 
            to : RegionSymbols
        ) -> dict:
    regionCoefficients = selectRegionCoefficients(from_, to)
    identifiers = [
            symbolicToIdentifier(regionCoefficients.transmissionCoefficent), 
            symbolicToIdentifier(regionCoefficients.reflectionCoefficient)
        ]
    identity = lambda identifier : sp.lambdify(identifiers, sp.Symbol(identifier))
    return {
            "from" : from_, 
            "to" : to, 
            "inputs" : identifiers, 
            "outputs" : identifiers, 
            "computations" : {
                    identifiers[0] : identity(identifiers[0]), 
                    identifiers[1] : identity(identifiers[1])
                }
        }

def makeCoefficentsFromHarmonicConstants(
            from_ : RegionSymbols, 
            to : RegionSymbols
        ) -> dict:
    regionCoefficients = selectRegionCoefficients(from_, to)
    harmonicConstants = [
            symbolicToIdentifier(from_.harmonicConstant), 
            symbolicToidentifier(to.harmonicConstant)
        ]
    transmission = 2 * from_.harmonicConstant \
            / (from_.harmonicConstant + to.harmonicConstant)
    reflection =  (from_.harmonicConstant - to.harmonicConstant) \
            / (from_.harmonicConstant + to.harmonicConstant)
    identifiers = [
            symbolicToIdentifier(regionCoefficients.transmissionCoefficent), 
            symbolicToIdentifier(regionCoefficients.reflectionCoefficient)
        ]
    return {
            "from" : from_, 
            "to" : to, 
            "inputs" : harmonicConstants, 
            "outputs" : identifiers, 
            "computations" : {
                    identifiers[0] : sp.lambdify(transmission, harmonicConstants), 
                    identifiers[1] : sp.lambdify(reflection, harmonicConstants)
                }
        }

def generateBoundryTransmissionReflectionCoefficients(
            regionSymbols : list[RegionSymbols], 
            transmissionReflectionCoefficentGenerator 
                    = defaultTransmissionReflectionCoefficientGenerator
        ) -> list: 
    return [
            transmissionReflectionCoefficentGenerator(
                    regionSymbols[ii + 1], 
                    regionSymbols[ii]
                ) \
            for ii in range(1, len(regionSymbols) - 1)
        ]

def calculateBoundryTransmissionReflectionCoefficients(
            calculators : list[dict], 
            inputs : dict
        ) -> dict: 
    results = []
    filterInputs = lambda calculator : {input_ : inputs[input_] for input_ in calculator["inputs"]}
    for calculator in calculators: 
        results.append({
                output : calculator["computations"][output](**filterInputs(calculator)) \
                for output in calculator["outputs"]
            })
    return results

def generateGeneralTransferFunctions(numberOfRegions : int) -> dict: 
    transferFunctionData = {"regionSymbols" : [RegionSymbols(ii - 1) for ii in range(numberOfRegions + 2)]}
    regionSymbols = transferFunctionData["regionSymbols"]
    transferFunctionData["transfers"] = [
            generalTransfer(
                    regionSymbols[ii + 1], 
                    regionSymbols[ii]
                ) \
            for ii in range(1, len(regionSymbols) - 1)
        ]
    transferFunctionData["transferFunctions"] = [
            lambdifyTransfer(transfer) for transfer in transferFunctionData["transfers"]
        ]
    return transferFunctionData

def constantPotentialTimeIndependentSchroedingerEquation1D( 
            regionSymbols : RegionSymbols, 
            totalEnergy : sp.Symbol = TOTAL_ENERGY_SYMBOL, 
            reducedPlanckConstant : sq.Quantity = hbar, 
            position : sp.Symbol = POSITION_SYMBOL 
        ): 
    return sp.Eq( 
            ((-(reducedPlanckConstant ** 2) / (2 * regionSymbols.mass)) \
                    * sp.Derivative(regionSymbols.waveFunction(position), (position, 2)))
                    + (regionSymbols.constantPotential * regionSymbols.waveFunction(position)), 
            totalEnergy * regionSymbols.waveFunction(position)
        )

def simpleWaveFunctionNormalization( 
            from_ : RegionSymbols, 
            toOrIndefinite : bool, 
            regionSymbols : RegionSymbols, 
            position : sp.Symbol = POSITION_SYMBOL, 
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

def manipulate(equation : sp.Eq, operation) -> sp.Eq: 
    return sp.Eq(operation(equation.lhs), operation(equation.rhs))

def secondDerivative(regionSymbols : RegionSymbols, position : sp.Symbol = POSITION_SYMBOL) -> sp.Derivative: 
    return sp.Derivative(regionSymbols.waveFunction(position), (position, 2))

def constantFactor(
            regionSymbols : RegionSymbols, 
            mass : sp.Symbol = MASS_SYMBOL, 
            reducedPlanckConstant : sq.Quantity = hbar
        ) -> sp.Mul:
    return (reducedPlanckConstant ** 2) / (2 * mass)

def secondDerivativeTerm(
            regionSymbols : RegionSymbols, 
            position : sp.Symbol = POSITION_SYMBOL, 
            mass : sp.Symbol = MASS_SYMBOL, 
            reducedPlanckConstant : sq.Quantity = hbar
        ) -> sp.Mul:
    return constantFactor(regionSymbols) * secondDerivative(regionSymbols)

def extractHarmonicConstant(regionSymbols : RegionSymbols, equation : sp.Eq): 
    waveFunction = regionSymbols.waveFunction(POSITION_SYMBOL)
    solveForWaveFunction = sp.Eq(waveFunction, sp.solve(equation, regionSymbols.waveFunction(POSITION_SYMBOL))[0])
    return sp.Eq(
            regionSymbols.harmonicConstant, 
            manipulate(solveForWaveFunction, \
                    lambda step : sp.sqrt((step / secondDerivative(regionSymbols)) ** (-1))).refine().rhs
        )


def extractAmplitudeCoefficents(exponentialEquation : sp.Eq, position : sp.Symbol = POSITION_SYMBOL): 
    C0 = sp.Wild("C0")
    C1 = sp.Wild("C1")
    harmonic = sp.Wild("k")
    results = exponentialEquation.match(C0 * sp.exp(harmonic * position) + C1 * sp.exp(-harmonic * position))
    
    return {
            "transmission" : results[C0], 
            "reflection" : results[C1], 
            "harmonicConstant" : results[harmonic]
        }

def createGeneralSolution(regionSymbols : RegionSymbols, waveEquation : sp.Eq, normalization : sp.Eq): 
    boundries = {
        regionSymbols.waveFunction(regionSymbols.startDistance) : regionSymbols.boundry, 
        regionSymbols.waveFunction(regionSymbols.distance) : regionSymbols.nextBoundry
    }
    harmonicConstant = extractHarmonicConstant(regionSymbols, waveEquation)
    generalSolution = sp.dsolve(waveEquation, ics = boundries)
    generalSolution = generalSolution.subs({harmonicConstant.rhs : harmonicConstant.lhs})
    return {
            "boundries" : boundries, 
            "harmonicConstantEquation" : harmonicConstant, 
            "generalSolution" : generalSolution
        }

def solveUnconstrainedParticularSolution(regionSymbols : RegionSymbols, generalSolution : sp.Eq) -> dict: 
    amplitudes = extractAmplitudeCoefficents(generalSolution["generalSolution"].rhs, POSITION_SYMBOL)
    exponential = generalSolution["generalSolution"].subs({
            amplitudes["transmission"] : regionSymbols.constants[0], 
            amplitudes["reflection"] : regionSymbols.constants[1]
        })
    return {
            "generalSolution" : generalSolution, 
            "exponential" : exponential, 
            "amplitudes" : amplitudes
        }

def extractNonZero(symbolic):
    harmonic = sp.Wild('k')
    expression = sp.Wild('C')
    nonZero = sp.Ne(harmonic, 0)
    otherExpression = sp.Wild('A')
    otherCondition = sp.Wild('B')
    otherStuff = sp.Wild('Q')
    extracted = symbolic.match(otherStuff + sp.Piecewise((expression, nonZero), (otherExpression, otherCondition)))
    return {
            "nonZero" : extracted[expression], 
            "otherwise" : extracted[otherExpression], 
            "rest" : extracted[otherStuff]
        }

def makeAmplitudeEquation(
            amplitudeConstant : sp.Symbol, 
            amplitudeSolutions : list, 
            refine : bool = False
        ) -> sp.Eq: 
    amplitudeExpression = (amplitudeSolutions[0] * amplitudeSolutions[1]).simplify()
    return sp.Eq(
            amplitudeConstant * amplitudeConstant, 
            amplitudeExpression.refine() if refine == True else amplitudeExpression
        )

def makeAmplitudeIntermediateEquations(
            unconstrainedParticularSolution : dict, 
            nonZeroNormalizationCase : sp.Eq, 
            regionSymbols : RegionSymbols
        ) -> dict:
    c0Solutions = sp.solve(nonZeroNormalizationCase,  regionSymbols.constants[0])
    c1Solutions = sp.solve(nonZeroNormalizationCase, regionSymbols.constants[1])
    intermediateSolutions = {}
    intermediateSolutions["unconstrainedAmplitudes"] = {
            "transmission" : c0Solutions, 
            "reflection" : c1Solutions
        }
    c0Equation = makeAmplitudeEquation(
            unconstrainedParticularSolution["amplitudes"]["transmission"], 
            c0Solutions
        )
    c1Equation = makeAmplitudeEquation(
            unconstrainedParticularSolution["amplitudes"]["reflection"], 
            c1Solutions
        )
    intermediateSolutions["intermediateSolutions"] = {
            "transmission" : c0Equation, 
            "reflection" : c1Equation
        }
    return intermediateSolutions

def solutionsEqualTo(equationToSolve : sp.Eq, toSolveFor : sp.Symbol) -> list: 
    return [sp.Eq(toSolveFor, solution) for solution in sp.solve(equationToSolve, toSolveFor)]

def makeConstrainedSolution(intermediateSolutions : dict, regionSymbols : RegionSymbols) -> dict:
    constrainedParticularSolution = {}
    constrainedParticularSolution["amplitudeConstrainedSolutions"] = {
            "transmission" : sp.solve(
                        intermediateSolutions["intermediateSolutions"]["reflection"], 
                        regionSymbols.constants[0]
                    ), 
            "reflection" : sp.solve(
                        intermediateSolutions["intermediateSolutions"]["transmission"], 
                        regionSymbols.constants[1]
                    )
        }
    constrainedParticularSolution["amplitudeConstrainedSolution"] = {
            "transmission" : makeAmplitudeEquation(
                    regionSymbols.constants[0], 
                    constrainedParticularSolution["amplitudeConstrainedSolutions"]["transmission"], 
                    True
                ), 
            "reflection" : makeAmplitudeEquation(
                    regionSymbols.constants[1], 
                    constrainedParticularSolution["amplitudeConstrainedSolutions"]["reflection"], 
                    True
                )
        }
    return constrainedParticularSolution

def solveForConstrainedAmplitudeCoefficients(
            regionSymbols : RegionSymbols, 
            normalization : sp.Eq, 
            unconstrainedParticularSolution : dict
        ) -> dict: 
    integration = normalization.subs({
            regionSymbols.waveFunction(POSITION_SYMBOL) \
                    : unconstrainedParticularSolution["exponential"].rhs
        }).doit()
    extractedConditions = extractNonZero(integration.lhs)
    nonZeroNormalizationCase = sp.Eq(
            extractedConditions["rest"] + extractedConditions["nonZero"], 
            integration.rhs
        )
    intermediateSolutions = makeAmplitudeIntermediateEquations(
            unconstrainedParticularSolution, 
            nonZeroNormalizationCase, 
            regionSymbols
        )
    constrainedParticularSolution = makeConstrainedSolution(intermediateSolutions, regionSymbols)
    unconstrainedParticularSolution["expandedExponential"] \
            = (unconstrainedParticularSolution["exponential"].rhs \
                    * sp.conjugate(unconstrainedParticularSolution["exponential"].rhs)
              ).expand()
    return {
            "regionSymbols" : regionSymbols, 
            "unconstrainedParticularSolution" : unconstrainedParticularSolution, 
            "intermediateSolutions" : intermediateSolutions, 
            "constrainedParticularSolution" : constrainedParticularSolution
        }

def makeBoundrySubstitutions(regionSymbols : RegionSymbols) -> dict: 
    boundrySubstitutions = {
            regionSymbols.boundry 
                    : symbolicToIdentifier(regionSymbols.waveFunction(regionSymbols.startDistance)), 
            regionSymbols.nextBoundry 
                    : symbolicToIdentifier(regionSymbols.waveFunction(regionSymbols.distance))
        }
    return boundrySubstitutions

def makeRootEquation(radicand, inverseCoefficent, scalar, offset): 
    return ((offset + sp.sqrt(radicand)) * scalar) / inverseCoefficent 

def extrapolateAmplitudeCoefficientComponents(
            constants : dict, 
            constrainedAmplitudeParticularSolutions : dict, 
            waveKeys : list = ["transmission", "reflection"], 
            radicandTag : str = RADICAND_TAG, 
            inverseCoefficentTag : str = INVERSE_COEFFICIENT_TAG, 
            offsetTag : str = OFFSET_TAG, 
            scalarTag : str = SCALAR_TAG
        ) -> dict: 
    radicand = sp.Wild(radicandTag)
    inverseCoefficent = sp.Wild(inverseCoefficentTag)
    offset = sp.Wild(offsetTag)
    scalar = sp.Wild(scalarTag)
    root = makeRootEquation(radicand, inverseCoefficent, scalar, offset)
    extrapolated = {}
    for wave in waveKeys: 
        constantName = symbolicToIdentifier(constants[wave])
        current = extrapolated[constantName] = {}
        current["components"] = []
        for solutionIndex, solution in enumerate(constrainedAmplitudeParticularSolutions[wave]): 
            rootComponents = solution.match(root)
            nameComponent = lambda component : constantName + str(solutionIndex) + str(component)[:-1]
            indexComponent = lambda componentName, component : {
                    "component" : sp.Abs(component) if componentName == radicandTag else component, 
                    "canonicalName" : str(componentName)[:-1]
                }
            rootComponents = {nameComponent(key) : indexComponent(key, item) for key, item in rootComponents.items()}
            current["components"].append(rootComponents)
    return extrapolated

def generateAmplitudCoefficientEquations(
            regionSymbols : RegionSymbols, 
            waveEquation : sp.Eq, 
            normalization : sp.Eq
        ) -> dict: 
    generalSolution = createGeneralSolution(regionSymbols, waveEquation, normalization)
    unconstrainedParticularSolution = solveUnconstrainedParticularSolution(regionSymbols, generalSolution)
    constrainedSolutions = solveForConstrainedAmplitudeCoefficients(
            regionSymbols, 
            normalization, 
            unconstrainedParticularSolution
        )
    constants = {
            "transmission" : regionSymbols.constants[0], 
            "reflection" : regionSymbols.constants[1]
        }
    extrapolatedParts = extrapolateAmplitudeCoefficientComponents(
            constants, 
            constrainedSolutions["constrainedParticularSolution"]["amplitudeConstrainedSolutions"]
        )
    return {
            "regionSymbols" : regionSymbols, 
            "waveEquation" : waveEquation, 
            "normalization" : normalization, 
            "generalSolution" : generalSolution, 
            "unconstrainedParticularSolution" : unconstrainedParticularSolution, 
            "constarinedParticularSolutions" : constrainedSolutions, 
            "extrapolatedComponentsOfConstants" : extrapolatedParts
        }

def lambdifySolutionComponents(solution : dict) -> dict: 
    componentFunctions = {}
    canoticalParameters = {}
    allParameters = set()
    for componentName, component in solution.items(): 
        componentExpression = component["component"]
        parameters = orderNames(list(
                substituteIdentifierAtomsList(componentExpression).values()
            ))
        allParameters |= set(parameters)
        name = functionNameFromIdentifier(componentName)
        componentFunctions[name] = {
                "function" : sp.lambdify(
                        parameters, 
                        componentExpression
                    ), 
                "parameters" : parameters
            }
        canoticalParameters[component["canonicalName"]] = name
    return {
            "componentFunctions" : componentFunctions, 
            "canoticalParameters" : canoticalParameters, 
            "allParameters" : allParameters
        }


def generateNumericalFunctionFromGeneralFunctionAndComponents(
            solution, 
            components : dict, 
            allParameters : dict
        ) -> dict: 
    def calculate(**kwargs): 
        generalFunctionArguments = {}
        for componentName, component in components.items(): 
            componentArguments = {parameter : kwargs[parameter] \
                    for parameter in component["parameters"]}
            generalFunctionArguments[componentName] \
                    = component["function"](**componentArguments)
        return solution["generalFunction"](**generalFunctionArguments)
    return calculate
    
def generateAmplitudeCoefficientNumericalFunctionsFromComponentEquations(
            amplitudeCoefficientEquations : dict, 
            radicandTag : str = RADICAND_TAG, 
            inverseCoefficentTag : str = INVERSE_COEFFICIENT_TAG, 
            offsetTag : str = OFFSET_TAG, 
            scalarTag : str = SCALAR_TAG
        ) -> dict: 
    constantComponents = amplitudeCoefficientEquations["extrapolatedComponentsOfConstants"]
    constantAmplitudeCalculations = {}
    for constant, components in constantComponents.items(): 
        solutions = []
        for solution in components["components"]: 
            numericalSolutionComponents = lambdifySolutionComponents(solution)
            allParameters = numericalSolutionComponents["allParameters"]
            solution["componentFunctions"] = numericalSolutionComponents["componentFunctions"]
            solution["generalFunctionParameters"] = orderNames(list(
                    numericalSolutionComponents["canoticalParameters"].values()
                ))
            canonicalParameters = numericalSolutionComponents["canoticalParameters"]
            canonicalParameters = {key : sp.Symbol(value) for key, value in canonicalParameters.items()}
            solution["generalFunction"] = sp.lambdify(
                    solution["generalFunctionParameters"], 
                    makeRootEquation(**canonicalParameters)
                )
            solution["parameters"] = orderNames(list(allParameters))
            solution["function"] = generateNumericalFunctionFromGeneralFunctionAndComponents(
                    solution, 
                    solution["componentFunctions"], 
                    allParameters
                )
            solution["numericalSolutionComponents"] = numericalSolutionComponents
            solutions.append(solution)
        constantAmplitudeCalculations[symbolicToIdentifier(constant)] = solutions
    return constantAmplitudeCalculations

def generateAmplitudeCoefficientNumericalFunctions(
            regionSymbols : RegionSymbols, 
            waveEquation : sp.Eq, 
            normalization : sp.Eq
        ) -> dict:
    amplitudeCoefficientEquations = generateAmplitudCoefficientEquations(
            regionSymbols, 
            waveEquation, 
            normalization
        )
    constantAmplitudeCalculations \
            = generateAmplitudeCoefficientNumericalFunctionsFromComponentEquations(
                    amplitudeCoefficientEquations
                )
    return {
            "amplitudeCoefficientEquations" : amplitudeCoefficientEquations, 
            "constantAmplitudeCalculations" : constantAmplitudeCalculations
        }

#def makeWaveFunctionFromTransfer(transferData : dict, transfers : dict) -> dict: 
#    calculation = transferData["transfer"].rhs
#    transferData["transfer"] 
#    nonNormalized = regionSymbols.constants[0] + regionSymbols.constants[1]
#    symbolicBoundry = (nonNormalized * sp.conjugate(nonNormalized)).simplify().refine()
#    outputs = transferData["outputs"]
#    nonNormalizedNumericWave = transfers[outputs[0]] + transfers[outputs[1]]
#    nonNormalizedNumericWave * np.conjugate(nonNormalizedNumericWave)
#    return {
#            "symbolic" : symbolicBoundry, 
#            "numeric" : nonNormalizedNumericWave
#        }


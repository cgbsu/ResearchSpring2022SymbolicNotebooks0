import types
from enum import Enum
from functools import partial
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal
from scipy import sparse
from scipy.sparse.linalg import eigsh

class EigenValueTypes(Enum): 
    LARGEST_MAGNITUDE = "LM"
    SMALLEST_MAGNITUDE = "SM"
    LARGEST_ALGEBRAIC = "LA"
    SMALLEST_ALGEBRAIC = "SA"
    HALF_SPECTRUM = "BE"

defaultScalingFactor : float = 1
defaultReducedPlancksConstant : float = 1
defaultMass : float = 1.0
defaultLength : float = 1.0
defaultUnormalizedPositionStep : float = 1e-3

def unormalizedPotentialTerm(
            potentialsOrPotentialFunction, 
            normalizedPositions : np.array, 
            length : float = defaultLength, 
            mass : float = defaultMass, 
            scalingFactor : float = defaultScalingFactor
        ) -> np.array: 
    if callable(type(potentialsOrPotentialFunction)) == True \
            and not hasattr(potentialsOrPotentialFunction, "__len__"): 
        potentials = potentialsOrPotentialFunction(normalizedPositions - length / 2.0)
    else:
        potentials : np.array = potentialsOrPotentialFunction
    return (length ** 2) * mass * potentials * scalingFactor

#def gridToDiagonal(grid : np.ndarray) -> np.ndarray: 
#    grid : np.array = grid.ravel()
#    return grid * np.identity(len(grid))

def positionGrid(pointCount : int) -> list[np.ndarray, np.ndarray]:
    return np.meshgrid(
            np.linspace(0, 1, pointCount, dtype=float), 
            np.linspace(0, 1, pointCount, dtype=float)
        )

def secondDerivativeOperator(
            pointCount : int, 
            reducedPlanckConstant : float = defaultReducedPlancksConstant, 
            mass : float = defaultMass
        ) -> np.ndarray: 
    ones : np.array = np.ones([pointCount])
    oneDimensionalDerivativeOperator = sparse.spdiags(
            np.array([ones, -2 * ones, ones]), 
            np.array([-1, 0, 1]), 
            pointCount, 
            pointCount
        )
    kineticEnegry = (-reducedPlanckConstant / (2.0 * mass)) * sparse.kronsum(
                oneDimensionalDerivativeOperator, 
                oneDimensionalDerivativeOperator
            )
    return \
        oneDimensionalDerivativeOperator, \
        kineticEnegry 

def makeHamiltonian(
            pointCount : int, 
            potential : np.ndarray, 
            reducedPlanckConstant : float = defaultReducedPlancksConstant, 
            mass : float = defaultMass
        ) -> np.ndarray: 
    oneDimensionalDerivativeOperator, kineticEnergy = secondDerivativeOperator(pointCount)
    return kineticEnergy + potential 

def compute2dWaveFunction(
            pointCount : int, 
            potential : np.ndarray, 
            eigenValueType : EigenValueTypes = EigenValueTypes.SMALLEST_MAGNITUDE, 
            energyCount : int = 10, 
            reducedPlanckConstant : float = defaultReducedPlancksConstant, 
            mass : float = defaultMass
        ) -> tuple[np.array, np.ndarray]: 
    reshapedPotential = sparse.diags(potential.reshape(pointCount ** 2), (0))
    hamiltonian : np.ndarray = makeHamiltonian(pointCount, reshapedPotential)
    energies, waveFunctions = eigsh(hamiltonian, k = energyCount, which = eigenValueType.value)
    waveFunctionGrids = np.array(list(map(
            lambda transposedWaveFunction : transposedWaveFunction.reshape((pointCount, pointCount)), 
            waveFunctions.T
        )))
    return energies, waveFunctionGrids 

def infiniteSquareWell(
            xPositions : np.ndarray, 
            yPositions : np.ndarray
        ) -> np.ndarray: 
    return 0 * xPositions

def gaussian(
            positions : np.ndarray, 
            offset : float, 
            variance : float
        ) -> np.ndarray: 
    exponent = (-1.0 / 2.0) * ((positions - offset) ** 2) / (2.0 * variance ** 2)
    scalar = 1.0 / (np.sqrt(2.0 * np.pi) * variance)
    return scalar * np.exp(exponent)
    

def gaussian2d(
            xPositions : np.ndarray, 
            yPositions : np.ndarray, 
            xOffset : float, 
            yOffset : float, 
            variance : float
        ) -> np.ndarray: 
    return gaussian(xPositions, xOffset, variance) * gaussian(yPositions, yOffset, variance)

def hydrogenAtom(
            xPosition : float, 
            yPosition : float, 
            centerX : float, 
            centerY : float, 
            bottom : float
        ) -> np.ndarray: 
    distance = 1 / (np.sqrt(((xPosition - centerX) ** 2) + ((yPosition - centerY) ** 2)) + bottom)
    return distance

def tunnelingCase(
            xPositions : np.ndarray, 
            yPositions : np.ndarray, 
            barrierPosition : float, 
            barrierWidth : float, 
            potentialHeight : float
        ) -> np.ndarray: 
    potentials = np.zeros(xPositions.shape)
    return np.where(
            (xPositions <= (barrierPosition + barrierWidth)) \
                    & (xPositions >= (barrierPosition - barrierWidth)), 
            potentialHeight, 
            potentials
        )

def finiteSquareWell(
            xPositions : np.ndarray, 
            yPositions : np.ndarray, 
            potentialHeight : float, 
            width : float, 
            height : float
        ) -> np.ndarray: 
    potential : np.ndarray = np.zeros(xPositions.shape)
    potential = np.where((xPositions <= width) | (xPositions >= (1 - width)), potentialHeight, potential)
    potential = np.where((yPositions <= height) | (yPositions >= (1 - height)), potentialHeight, potential)
    return potential

def finiteCircularWell(
            xPositions : np.ndarray, 
            yPositions : np.ndarray, 
            potentialHeight : float, 
            radius : float, 
            centerX : float = .5, 
            centerY : float = .5
        ) -> np.ndarray: 
    distances : float = np.sqrt(((xPositions - centerX) ** 2) + ((yPositions - centerY) ** 2))
    potential : np.ndarray = np.zeros(xPositions.shape)
    potential = np.where(distances >= radius, potentialHeight, potential)
    return potential

def stairwell(
            xPositions : np.ndarray, 
            yPositions : np.ndarray, 
            unitPotentialHeight : float, 
            widthRatios : list[float], 
            heightRatios : list[float], 
            unitLength : float  = 1
        ) -> np.ndarray: 
    assert len(heightRatios) == len(widthRatios)
    potential : np.ndarray = np.zeros(xPositions.shape)
    previousLength : float = 0
    lengthRatio : float = 0
    for ii in range(len(widthRatios)): 
        lengthRatio += widthRatios[ii]
        length = lengthRatio * unitLength
        potential = np.where(
                (xPositions <= (length))
                        & (xPositions >= previousLength), 
                heightRatios[ii] * unitPotentialHeight, 
                potential
            )
        previousLength = length
    return potential


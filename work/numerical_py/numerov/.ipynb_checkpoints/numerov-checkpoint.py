import types

from functools import partial

import sympy as sp

import numpy as np

import matplotlib.pyplot as plt

from scipy.linalg import eigh_tridiagonal



# Following: https://youtu.be/ay0zZ8SUMSk

# Find psi_1 to psi_{N - 1}, excluding psi_0 = 0 and psi_N = 0



defaultScalingFactor : float = 1000

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

def makeEigenMatrixTerms(
            potential : np.array, 
            normalizedStep : float
        ) -> tuple[np.array, np.array]: 
    inverseNormalizedStepSquared : float = 1.0 / (normalizedStep ** 2)
    return (potential + inverseNormalizedStepSquared, -inverseNormalizedStepSquared / 2.0 * np.ones(len(potential)))

def makeEigenFunctions(
            potential : np.array, 
            normalizedStep : float, 
        ) -> tuple[np.array, np.ndarray]: 
    pointCount = len(potential)
    assert pointCount > 2, "Need at least 3 points to create eigen matrix"
    potentialTerm, stepTerm = makeEigenMatrixTerms(potential, normalizedStep)
    return eigh_tridiagonal(potentialTerm[ : -1], stepTerm[1 : -1])



def finiteSquareWell(
            normalizedPosition : np.array, 
            potentialHeight : float, 
            length : float
        ) -> np.array: 
    potentials = np.zeros(len(normalizedPosition))
    potentials =  np.where(
            normalizedPosition > (length / 3.0), 
            potentials, 
            potentialHeight
        )
    potentials =  np.where(
            normalizedPosition < (2.0 * length / 3.0), 
            potentials, 
            potentialHeight
        )
    return potentials

def tunnelingCase(
            normalizedPosition : np.array, 
            potentialHeight : float, 
            length : float
        ) -> np.array: 
    potentials = np.zeros(len(normalizedPosition))
    potentials =  np.where(
            ~((normalizedPosition > (2.0  * length / 4.0)) \
                & (normalizedPosition < (3.0  * length / 4.0))), 
            potentials, 
            potentialHeight
        )
    return potentials

def computeWaveFunction(
            potentialsOrPotentialFunction, 
            scalingFactor : float = defaultScalingFactor, 
            positionStep : float = defaultUnormalizedPositionStep, 
            length : float = defaultLength, 
            mass : float = defaultMass, 
            reducedPlanckConstant : float = defaultReducedPlancksConstant
        ) -> dict[np.array, np.array, np.ndarray]:
    normalizedPositions : np.array = np.arange(0, length, positionStep) / length  
    normalizedPositionStep = 1.0 / len(normalizedPositions)
    potential = unormalizedPotentialTerm(
            potentialsOrPotentialFunction, 
            normalizedPositions, 
            length, 
            mass, 
            scalingFactor
        )
    energies, waveFunctions = makeEigenFunctions(potential, normalizedPositionStep)
    return {
            "normalizedPositions" : normalizedPositions, 
            "potential" : potential, 
            "energies" : energies, 
            "waveFunctions" : waveFunctions.T
        }

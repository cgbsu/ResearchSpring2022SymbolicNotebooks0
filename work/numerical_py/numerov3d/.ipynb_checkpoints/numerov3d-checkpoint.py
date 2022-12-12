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
class DimensionIndex(Enum):
    X = 0
    Y = 1
    Z = 2
class WaveFunctions:
    def __init__(
                self, 
                shape : tuple[int], 
                energyValues : np.array, 
                eigenVectors : np.ndarray, 
                doNotComputeExtra : bool = False
            ):
        self.shape : tuple[int] = shape
        self.pointCount = shape[0]
        self.dimensions = len(shape)
        self.energyValues : np.array = energyValues
        self.waveFunctions : np.ndarray = np.array(list(map(
                lambda transposedWaveFunction : transposedWaveFunction.reshape(self.shape), 
                eigenVectors.T
            )))
        if doNotComputeExtra == True: 
            self.probabilities = None
            self.decibleProbabilities = None
        else: 
            self.probabilities = self.waveFunctions * np.conjugate(self.waveFunctions)
            self.decibleProbabilities = 10 * np.log10(self.probabiliti
class MeshGrid: 
    def __init__(self, gridDimensionalComponents : tuple[np.ndarray], pointCount : int, length : float): 
        self.pointCount = pointCount
        self.length = length
        self.gridDimensionalComponents : tuple[np.ndarray] = gridDimensionalComponents 
        self.dimensions = len(self.gridDimensionalComponents)
        for dimension_ in list(DimensionIndex.__members__): 
            dimension = getattr(DimensionIndex, dimension_)
            if self.dimensions > dimension.value: 
                setattr(self, dimension.name.lower(), self.gridDimensionalComponents[dimension.valu
def makeLinspaceGrid(pointCount : int, length : float, dimensions : int, componentType : type = float) -> MeshGrid: 
    spaces : tuple[np.array] = tuple((np.linspace(0, length, pointCount, dtype = componentType) for ii in range(dimensions)))
    return MeshGrid(np.meshgrid(*spaces), pointCount, leng
def makeMappingMatrix(pointCount : int, dimensions : int) -> np.ndarray:
    ones = np.ones([pointCount])
    baseMappingMatrix = sparse.spdiags(
        np.array([ones, -2 * ones, ones]), 
        np.array([-1, 0, 1]), 
        pointCount, 
        pointCount
    )
    mappingMatrix = baseMappingMatrix
    for ii in range(1, dimensions): 
        mappingMatrix = sparse.kronsum(mappingMatrix, baseMappingMatrix)
    return mappingMatrix, baseMappingMat
def kineticEnergyOperator(mappingMatrix : np.ndarray) -> np.ndarray: 
    return (-1.0 / 2.0) * mappingMat
def potentialEnergyOperator(potential : np.ndarray, pointCount : int, dimensions : int) -> np.ndarray: 
    return sparse.diags(potential.reshape(pointCount ** dimensions), (
def makeHamiltonian(
            potential : np.ndarray, 
            pointCount : int, 
            dimensions : int, 
            mappingMatrix : np.ndarray
        ) -> np.ndarray: 
    return kineticEnergyOperator(mappingMatrix) + potentialEnergyOperator(
            potential, 
            pointCount, 
            dimensions
def computeWaveFunction(
            potential : np.ndarray, 
            energyCount : int = 10, 
            eigenValueType : EigenValueTypes = EigenValueTypes.SMALLEST_MAGNITUDE
        ) -> WaveFunctions: 
    dimensions : int = len(potential.shape)
    pointCount : int = potential.shape[0]
    for cardinality in potential.shape: 
        assert cardinality == pointCount, "All dimensions of potential need to have the same number of elements"
    mappingMatrix, _ = makeMappingMatrix(pointCount, dimensio
    return WaveFunctions(potential.shape, *eigsh(
            makeHamiltonian(potential, pointCount, dimensions, mappingMatrix), 
            k = energyCount, 
            which = eigenValueType.value
        )


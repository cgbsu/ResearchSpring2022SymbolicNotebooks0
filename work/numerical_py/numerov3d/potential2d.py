import numpy as np

def gaussian(
            positions : np.ndarray, 
            offset : float, 
            variance : float
        ) -> np.ndarray: 
    exponent = (-1.0 / 2.0) * ((positions - offset) ** 2) / (2.0 * variance ** 2)
    scalar = 1.0 / (np.sqrt(2.0 * np.pi) * variance)
    return scalar * np.exp(exponent)
    

def gaussian2d(
            grid, 
            xOffset : float, 
            yOffset : float, 
            variance : float
        ) -> np.ndarray: 
    return gaussian(grid.x, xOffset, variance) * gaussian(grid.x, yOffset, variance)

def hydrogenAtom(
            grid, 
            centerX : float, 
            centerY : float, 
            bottom : float
        ) -> np.ndarray: 
    distance = 1 / (np.sqrt(((xPosition - centerX) ** 2) + ((yPosition - centerY) ** 2)) + bottom)
    return distance

def tunnelingCase(
            grid, 
            barrierPosition : float, 
            barrierWidth : float, 
            potentialHeight : float
        ) -> np.ndarray: 
    potentials = np.zeros(xPositions.shape)
    return np.where(
            (grid.x <= (barrierPosition + barrierWidth)) \
                    & (grid.x >= (barrierPosition - barrierWidth)), 
            potentialHeight, 
            potentials
        )

def finiteSquareWell(
            grid, 
            potentialHeight : float, 
            width : float, 
            height : float
        ) -> np.ndarray: 
    potential : np.ndarray = np.zeros(xPositions.shape)
    potential = np.where((grid.x <= width) | (grid.x >= (1 - width)), potentialHeight, potential)
    potential = np.where((grid.y <= height) | (grid.y >= (1 - height)), potentialHeight, potential)
    return potential

def finiteCircularWell(
            grid, 
            potentialHeight : float, 
            radius : float, 
            centerX : float = .5, 
            centerY : float = .5
        ) -> np.ndarray: 
    distances : float = np.sqrt(((grid.x - centerX) ** 2) + ((grid.y - centerY) ** 2))
    potential : np.ndarray = np.zeros(grid.x.shape)
    potential = np.where(distances >= radius, potentialHeight, potential)
    return potential

def stairwell(
            grid, 
            unitPotentialHeight : float, 
            widthRatios : list[float], 
            heightRatios : list[float], 
            unitLength : float  = 1
        ) -> np.ndarray: 
    assert len(heightRatios) == len(widthRatios)
    potential : np.ndarray = np.zeros(grid.x.shape)
    previousLength : float = 0
    lengthRatio : float = 0
    for ii in range(len(widthRatios)): 
        lengthRatio += widthRatios[ii]
        length = lengthRatio * unitLength
        potential = np.where(
                (grid.x <= (length))
                        & (grid.x >= previousLength), 
                heightRatios[ii] * unitPotentialHeight, 
                potential
            )
        previousLength = length
    return potential


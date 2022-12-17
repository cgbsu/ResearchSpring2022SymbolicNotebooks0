import gdspy
from picwriter import toolkit as tk
import picwriter.components as pc
from enum import Enum

class FeatureType(Enum): 
    ABOVE_SURFACE = "above_surface"
    ON_SURFACE = "ignore_on_surface"
    BELOW_SURFACE = "negative_space_below_suraface"

class RectangleGrid2DTemplate: 
    DEFAULT_LAYERS = {
           FeatureType.ABOVE_SURFACE : 102, 
           FeatureType.ON_SURFACE : 101, 
           FeatureType.BELOW_SURFACE : 100
        }
    def __init__(
                self, 
                lattice_constant, 
                radius, 
                layers = DEFAULT_LAYERS, 
                feature_type = FeatureType.ABOVE_SURFACE, 
                defect_type = FeatureType.ON_SURFACE, 
                ignore_defects = True, 
                datatype = 0
            ): 
       self.lattice_constant = lattice_constant
       self.radius = radius
       self.feature_type = feature_type
       self.defect_type = defect_type
       self.ignore_defects = ignore_defects
       self.layers = layers
       self.datatype = datatype


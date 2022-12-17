from enum import Enum
import gdspy
from picwriter import toolkit as tk
import picwriter.components as pc
import custom_components.RectangularGrid2DTemplate as rg2dt
import numpy as np

class DimensionalIndex(Enum): 
    X = 0
    Y = 1
    Z = 2

class RectangularGrid2D(tk.Component): 
    GRID_DEFECT_FEATURE_MAPPING = {
            "feature" : 1, 
            "defect" : -1
        }
    FULL_AXIS = -1
    GRID_DEFECT_FEATURE = GRID_DEFECT_FEATURE_MAPPING["defect"]
    def __init__(
                self, 
                template, 
                extent, 
                defect_list = [], 
                position = (0, 0), 
                resist = "+",
                fab = "ETCH"
            ): 
        tk.Component.__init__(self, "RectangularGrid2D", locals())
        self.template = template
        self.position = position
        self.extent = extent
        self.defect_list = defect_list
        self.resist = resist
        self.fab = fab

        self.__build_grid()
        self.__build_cell()
        self.__build_ports()

        """ Translate & rotate the ports corresponding to this specific component object
        """
        self._auto_transform_()

    def __build_grid(self): 
        assert RectangularGrid2D.GRID_DEFECT_FEATURE_MAPPING["feature"] == 1
        GRID_DEFECT_FEATURE = RectangularGrid2D.GRID_DEFECT_FEATURE 
        grid = np.ones((
                self.extent[DimensionalIndex.X.value], 
                self.extent[DimensionalIndex.Y.value]
            ))
        for defect_indicies in self.defect_list: 
            x_index = defect_indicies[DimensionalIndex.X.value]
            y_index = defect_indicies[DimensionalIndex.Y.value]
            if x_index == RectangularGrid2D.FULL_AXIS: 
                grid[:, y_index] = GRID_DEFECT_FEATURE
            elif y_index == RectangularGrid2D.FULL_AXIS:
                grid[x_index, :] = GRID_DEFECT_FEATURE
            else: 
                grid[x_index, y_index] = GRID_DEFECT_FEATURE
        self.grid = grid

    def __build_cell(self): 
        GRID_DEFECT_FEATURE = RectangularGrid2D.GRID_DEFECT_FEATURE 
        radius = self.template.lattice_constant * self.template.radius
        x_position = radius + self.position[DimensionalIndex.X.value]
        y_position = radius + self.position[DimensionalIndex.Y.value]
        print("Layers: ", self.template.layers, "\nFEATURE TYPE: ", self.template.feature_type)
        for x_index in range(self.extent[DimensionalIndex.X.value]): 
            x = self.template.lattice_constant * x_index + x_position
            for y_index in range(self.extent[DimensionalIndex.Y.value]): 
                y = self.template.lattice_constant * y_index + y_position
                if self.grid[x_index][y_index] != GRID_DEFECT_FEATURE: 
                    self.add(
                            gdspy.Round((x, y), radius), 
                            self.template.layers[self.template.feature_type], 
                            self.template.datatype
                        )
                elif self.template.ignore_defects == False: 
                    self.add(
                            gdspy.Round((x, y), radius), 
                            self.template.layers[self.template.defect_type], 
                            self.template.datatype
                        )
        print("Final position: ", x, y)

    def __build_ports(self): 
        self.portlist["output"] = {"port": (0, 0), "direction": "EAST"}

if __name__ == "__main__": 
    waveGuideTemplate = pc.WaveguideTemplate()
    waveguide = pc.Waveguide([(0, 500), (1000, 500)], waveGuideTemplate)
    gridTemplate = rg2dt.RectangleGrid2DTemplate(1, .3)
    plainGrid = RectangularGrid2D(gridTemplate, (100, 100))
    columnMissingGrid = RectangularGrid2D(gridTemplate, (100, 100), [(10, -1)], (0, 100))
    rowMissingGrid = RectangularGrid2D(gridTemplate, (100, 100), [(-1, 10)], (200, 100))
    sparseDefectGrid = RectangularGrid2D(gridTemplate, (100, 100), [(10, 10), (40, 32), (28, 34), (3, 45)], (300, 100))
    top = gdspy.Cell("top")
    tk.add(top, waveguide)
    tk.add(top, plainGrid)
    tk.add(top, columnMissingGrid)
    tk.add(top, rowMissingGrid)
    tk.add(top, sparseDefectGrid)
    tk.build_mask(top, waveGuideTemplate, final_layer = 200, final_datatype = 0)
    gdspy.LayoutViewer()



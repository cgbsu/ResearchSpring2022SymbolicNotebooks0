from enum import Enum
import sys
from pathlib import Path
import numpy as np
import picwriter.components as pc
from picwriter import toolkit as tk
sys.path.insert(0, '..')
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
import custom_components as cc

class Vector: 
    def __init__(self, vector_tuple : tuple[float, float]): 
        self.x = vector_tuple[cc.DimensionalIndex.X.value]
        self.y = vector_tuple[cc.DimensionalIndex.Y.value]
        self.dimensions = 2
        if len(vector_tuple) >= 3: 
            self.z = vector_tuple[cc.DimensionalIndex.Z.value]
            self.dimensions = 3

    def to_tuple(self) -> tuple: 
        if self.dimensions < 3:
            return (self.x, self.y)
        else: 
            return (self.x, self.y, self.z)

# Waveguide
# Input Pads, Output Pads
# Stairwell Pad Group(s)
# Wireing scheme

# Pad Group: 
#   - Pad (longitudonal) Position(s)
#   - Pad Input
#   - Pad Output
#   ~ Pad Template

class PadType(Enum): 
    TOP = 1, 
    BOTTOM = -1

class FixedStaggeredPadGroupTemplate: 
    def __init__(
                self, 
                template, 
                seperation : float, 
                width : float, 
                relative_input_locations : list[tuple[float, float]], 
                pad_to_wave_guide_offset : float, 
                potential_ratios = [1, 2 / 3, 1 / 3], 
                length_ratios = [1 / 3, 1 / 3, 1 / 3]
            ): 
        self.relative_input_locations = relative_input_locations
        self.template = template
        self.seperation = seperation, 
        self.width = width
        self.potential_ratios = potential_ratios
        self.length_ratios = length_ratios
        self.pad_to_wave_guide_offset = pad_to_wave_guide_offset

def simple_build_trace(
            from_port : tuple[float, float], 
            to_port : tuple[float, float]
        ) -> list[tuple[float, float]]: 
    from_port = Vector(from_port)
    to_port = Vector(to_port)
    if from_port.x == to_port.x or from_port.y == to_port.y: 
        return [from_port.to_tuple(), to_port.to_tuple()]
    else: 
        return [
                from_port.to_tuple(), 
                (from_port.x, to_port.y), 
                to_port.to_tuple()
            ]

class PadGroup: 
    def __init__(self, pad, routs): 
        self.pad = pad
        self.routs = routs
    def integrate(self, top): 
        tk.add(top, self.pad)
        for route in self.routs: 
            tk.add(top, route)

def total_stairwell_lifetime_scaled_pad_length(
            wave_length : float, 
            mean_life_time : float, 
            speed_of_light : float, 
            lifetime_to_bondpad_length_ratio : float = 3, 
            scale_to_wave_length = False
        ): 
    return mean_life_time * speed_of_light * lifetime_to_bondpad_length_ratio \
            if scale_to_wave_length == False \
            else wave_length * lifetime_to_bondpad_length_ratio

def place_lifetime_scaled_static_pad(
            wave_guide_center_y : float, 
            x : float, 
            total_length : float, 
            wave_guide_width : float, 
            template : FixedStaggeredPadGroupTemplate, 
            pad_type : PadType, 
            bond_pad_ports : list[tuple[float, float]], 
            metal_route_template : pc.MetalTemplate, 
            mean_life_time : float, 
            speed_of_light : float, 
            connection_matrix : tuple[tuple[bool]] = ((True, False), ), 
            lifetime_to_bondpad_length_ratio : float = 3, 
            scale_to_wave_length = False, 
            wave_length = 1.5, 
        ) -> PadGroup: 
    print("R: ", lifetime_to_bondpad_length_ratio)
    total_stairwell_pad_length = total_stairwell_lifetime_scaled_pad_length(
            wave_length, 
            mean_life_time, 
            speed_of_light, 
            lifetime_to_bondpad_length_ratio, 
            scale_to_wave_length
        )
    if (x + total_stairwell_pad_length) >= total_length: 
        assert total_length > total_stairwell_pad_length
        x = (total_length - total_stairwell_pad_length) / 2
    max_cladding_width = template.template.clad_width \
            * np.max(np.array(template.potential_ratios))
    y_offset = (template.width / 2) \
                    + template.pad_to_wave_guide_offset \
                    + max_cladding_width + wave_guide_width / 2
    offset_scalar = 1 if pad_type == PadType.TOP else -1
    pad = cc.StaggeredBondpad(
            template.template, 
            total_stairwell_pad_length, 
            template.width, 
            template.potential_ratios, 
            template.length_ratios, 
            port = (x, wave_guide_center_y + (offset_scalar * y_offset))
        )
    pad_position = pad.portlist["output"]["port"]
    routes = []
    for rout_input_index in range(len(template.relative_input_locations)): 
        for bond_pad_index in range(len(bond_pad_ports)): 
            if connection_matrix[bond_pad_index][rout_input_index] == True:
                routes.append(pc.MetalRoute(
                        simple_build_trace(
                                (
                                        pad.port[cc.DimensionalIndex.X.value] \
                                                + template.relative_input_locations[rout_input_index][cc.DimensionalIndex.X.value] \
                                                * total_stairwell_pad_length, 
                                        pad.port[cc.DimensionalIndex.Y.value] \
                                                + template.relative_input_locations[rout_input_index][cc.DimensionalIndex.Y.value] \
                                                * pad.maxWidth
                                ), 
                                bond_pad_ports[-1]
                            ), 
                        metal_route_template
                    ))
    return PadGroup(pad, routes)


def get_wave_guide_width(wave_guide) -> float: 
    if isinstance(wave_guide, pc.Waveguide): 
        return wave_guide.wgt.wg_width
    else: 
        extent = wave_guide.extent
        return (extent[1] - 1) * wave_guide.template.lattice_constant

def get_wave_guide_center_y_position(wave_guide) -> float: 
    if isinstance(wave_guide, pc.Waveguide): 
        return wave_guide.portlist["input"]["port"][cc.DimensionalIndex.Y.value]
    else: 
        return wave_guide.position[cc.DimensionalIndex.Y.value] + get_wave_guide_width(wave_guide) / 2

class Stairwell: 
    DEFAULT_START_X_RATIO = 1 / 5
    DEFAULT_GROUP_COUNT = 1
    DEFAULT_CONNECTION_MATRIX = {
            PadType.TOP : ((True, False), ), 
            PadType.BOTTOM : ((True, False), )
        }
    def __init__(
                self, 
                top, 
                wave_guide, 
                pad_group_template, 
                pad_group_builder, 
                total_length, 
                pad_group_count = DEFAULT_GROUP_COUNT , 
                start_x_ratio = DEFAULT_START_X_RATIO, 
                top_bond_pad_ports = [], 
                bottom_bond_pad_ports = [], 
                route_template = pc.MetalTemplate(), 
                builder_arguments = {}, 
                connection_matrix : tuple[tuple[bool]] = DEFAULT_CONNECTION_MATRIX
            ):
        self.top = top
        self.wave_guide = wave_guide
        self.pad_group_template = pad_group_template
        self.pad_group_builder = pad_group_builder
        self.wave_guide_width = get_wave_guide_width(self.wave_guide)
        self.wave_guide_center_y = get_wave_guide_center_y_position(self.wave_guide)
        self.pad_group_count = pad_group_count 
        self.builder_arguments = builder_arguments
        self.start_x_ratio = start_x_ratio
        self.total_length = total_length
        self.top_groups = []
        self.bottom_groups = []
        self.top_bond_pad_ports = top_bond_pad_ports 
        self.bottom_bond_pad_ports = bottom_bond_pad_ports 
        self.route_template = route_template
        self.portlist = {}
        self.connection_matrix = connection_matrix 

        self.__build_cell()
        self.__build_ports()

    def __build_cell(self): 
        x = self.start_x_ratio * self.total_length
        make_pad = lambda x, pad_type, bondpads : self.pad_group_builder(
                self.wave_guide_center_y, 
                x, 
                self.total_length, 
                self.wave_guide_width, 
                self.pad_group_template, 
                pad_type, 
                bondpads, 
                self.route_template, 
                **self.builder_arguments, 
                connection_matrix = self.connection_matrix[pad_type]
            )
        for ii in range(self.pad_group_count): 
            self.top_groups.append(make_pad(x, PadType.TOP, self.top_bond_pad_ports))
            print("T: ", self.top_bond_pad_ports)
            self.bottom_groups.append(make_pad(x, PadType.BOTTOM, self.bottom_bond_pad_ports))
            print("B: ", self.bottom_bond_pad_ports)
            x = self.top_groups[-1].pad.padExtents[-1]
            self.top_groups[-1].integrate(self.top)
            self.bottom_groups[-1].integrate(self.top)

    def __build_ports(self): 
        self.portlist["output"] = {"port": (0, 0), "direction": "EAST"}

    def vertical_width(self): 
        vertical_extents = []
        for group in self.top_groups[-1]: 
            vertical_extents.append(
                    group.pad.padExtents[-1][cc.DimensionalIndex.Y.value]
                )
        for group in self.bottom_groups[-1]: 
            vertical_extents.append(
                    group.pad.padExtents[-1][cc.DimensionalIndex.Y.value]
                )
        vertical_extents = np.array(vertical_extents)
        return np.abs(np.min(vertical_extents) - np.max(vertical_extents))

if __name__ == "__main__": 
    import gdspy
    OBLIGITORY_WAVE_GUIDE_TEMPLATE = pc.WaveguideTemplate()
    OBLIGITORY_WAVE_GUIDE = pc.Waveguide([(0, 0), (1, 0)], OBLIGITORY_WAVE_GUIDE_TEMPLATE)
    metal_template = pc.MetalTemplate()
    staggered_metal_template = cc.StaggeredMetalTemplate()
    total_length = 10000
    lattice_constant : float = .5
    radius : float = .2 * lattice_constant
    extent_y_index : int = 7
    extent_x_index = int(round(total_length / lattice_constant))
    defect_layer_index = int(np.floor(extent_y_index / 2))
    grid_template = cc.RectangleGrid2DTemplate(lattice_constant, radius)
    wave_guide_grid = cc.RectangularGrid2D(
            grid_template, 
            (extent_x_index, extent_y_index), 
            [(-1, defect_layer_index)], 
            position = (0, 500)
        )
    bottom_pad = pc.Bondpad(metal_template)
    top_pad = pc.Bondpad(metal_template, port=(0, 1000))
    top = gdspy.Cell("top")
    tk.add(top, wave_guide_grid)
    tk.add(top, top_pad)
    tk.add(top, bottom_pad)
    stairwell = Stairwell(
            top, 
            wave_guide_grid, 
            FixedStaggeredPadGroupTemplate(
                    staggered_metal_template, 
                    0, 
                    100, 
                    [], 
                    10
                ), 
            place_lifetime_scaled_static_pad, 
            total_length, 
            Stairwell.DEFAULT_GROUP_COUNT, 
            Stairwell.DEFAULT_START_X_RATIO, 
            [top_pad.portlist["input"]["port"]], 
            [bottom_pad.portlist["input"]["port"]], 
            metal_template, 
            builder_arguments = {
                    "speed_of_light" : 3e14, 
                    "mean_life_time" : 10e-12, 
                    "lifetime_to_bondpad_length_ratio" : 3
                }
        )
    tk.add(top, OBLIGITORY_WAVE_GUIDE)
    tk.build_mask(top, OBLIGITORY_WAVE_GUIDE_TEMPLATE, final_layer = 200, final_datatype = 0)
    gdspy.LayoutViewer()


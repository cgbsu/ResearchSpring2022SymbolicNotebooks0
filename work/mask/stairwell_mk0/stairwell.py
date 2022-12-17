from enum import Enum
import sys
from pathlib import Path
import numpy as np
import gdspy
from picwriter import toolkit as tk
import picwriter.components as pc

sys.path.insert(0, '..')
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
import custom_components as cc
#from custom_components.RectangularGrid2DTemplate import RectangleGrid2DTemplate#RectangularGrid2DTemplate#  as rg2dt
#from custom_components import RectangularGrid2D# as rg2d

MICROMETERS = 1e-6
UNITS = MICROMETERS # Micrometers
NANOMETERS_IN_METERS = 1e-9
SPEED_OF_LIGHT_METERS_PER_SECOND = 299792458
SPEED_OF_LIGHT_MICRO_METERS_PER_SECOND = SPEED_OF_LIGHT_METERS_PER_SECOND / UNITS
MICROMETERS_IN_MILLIMETER = 1000
CHIP_TOTAL_WIDTH = MICROMETERS_IN_MILLIMETER * 10 # .5 centimeters

OBLIGITORY_WAVE_GUIDE_TEMPLATE = pc.WaveguideTemplate()
OBLIGITORY_WAVE_GUIDE = pc.Waveguide([(0, 0), (1, 0)], OBLIGITORY_WAVE_GUIDE_TEMPLATE)

APPROXIMATLEY_LINEAR_FREQUENCY_REGION = (.35, .42) # Photonic Crystal's Modeling the Flow of Light, Chapter 10 Page 194
                                                   # Wikipedia "Group Velocity" https://en.wikipedia.org/wiki/Group_velocity
                                                   # Various papers doubting superluminal phenomonia/the Hartman Effect citing the change in the shape of the wave packet.

def calculate_lattice_constant_in_units(
            wavelength_in_units : float, 
            frequency : float, 
            units = UNITS
        ) -> float:
    # wavelength * lattice_constant / (2 * pi * c) = f
    # --> f * 2 * pi * c / wavelength = lattice_constatnt
    return frequency * wavelength_in_units

def solve_for_radius(lattice_constant : float, radius_over_lattice_constant_ratio : float) -> float: 
    return lattice_constant * radius_over_lattice_constant_ratio

def single_stairwell_rods(
            top, 
            lattice_constant : float, 
            radius : float, 
            mean_life_time : float, 
            lifetime_to_bondpad_length_ratio : float = 3, 
            speed_of_light = SPEED_OF_LIGHT_MICRO_METERS_PER_SECOND, 
            potential_height_ratios = [1, 2 / 3, 1 / 3], 
            segment_length_ratios = [1 / 3, 1 / 3, 1 / 3], 
            total_length = CHIP_TOTAL_WIDTH, 
            width = 7, 
            staggered_pad_width = 100, 
            offset_from_staggered_pad = 10
        ): 
    total_stairwell_pad_length = mean_life_time * speed_of_light * lifetime_to_bondpad_length_ratio
    print("Stairwell Pad Length: ", total_stairwell_pad_length)
    print("Total Length: ", total_length)
    staggered_pad_template = cc.StaggeredMetalTemplate()
    grid_template = cc.RectangleGrid2DTemplate(lattice_constant, radius)
    top_stairwell_pad_position = staggered_pad_width * 2
    left_stairwell_pad_position = total_length / (len(segment_length_ratios) + 2)
    if (left_stairwell_pad_position + total_stairwell_pad_length) > total_length: 
        assert total_stairwell_pad_length < total_length, "Stairwell too long for this chip (assuming speed of light in vaccume)!"
        left_stairwell_pad_position = (total_length - total_stairwell_pad_length) / 2
    extent_x_index = int(round(total_length / lattice_constant))
    defect_layer_index = int(np.floor(width / 2))
    print("Stairwell Left Pad Position: ", left_stairwell_pad_position)
    print("Total length: ", total_length)
    print("Extent Count: ", extent_x_index)
    print("Defect Layer Index: ", defect_layer_index)
    anode_pad = cc.StaggeredBondpad(
            staggered_pad_template, 
            total_stairwell_pad_length, 
            staggered_pad_width, 
            potential_height_ratios, 
            segment_length_ratios, 
            port = (left_stairwell_pad_position, top_stairwell_pad_position)
        )
    cathode_pad = cc.StaggeredBondpad(
            staggered_pad_template, 
            total_stairwell_pad_length, 
            staggered_pad_width, 
            potential_height_ratios, 
            segment_length_ratios, 
            port = (
                    left_stairwell_pad_position, 
                    top_stairwell_pad_position - staggered_pad_width \
                            - (anode_pad.maxCladdingWidth * 2) \
                            - (offset_from_staggered_pad * 2) - defect_layer_index
                )
        )
    wave_guide_grid = cc.RectangularGrid2D(
            grid_template, 
            (extent_x_index, width), 
            [(-1, defect_layer_index)], 
            position = (
                    0, 
                    top_stairwell_pad_position \
                            - (staggered_pad_width / 2) \
                            - anode_pad.maxCladdingWidth \
                            - defect_layer_index \
                            - offset_from_staggered_pad  
                )
        )
    tk.add(top, wave_guide_grid)
    tk.add(top, anode_pad)
    tk.add(top, cathode_pad)

def main(): 
    wave_length = 1550 * (NANOMETERS_IN_METERS / UNITS)
    radius_to_lattice_constant_ratio = .2 # TODO: Update from .2 for concodinanite. Photonic Crystal's Modeling the Flow of Light, Chapter 10 Page 192
    mean_life_time_estimation = 10e-12 # TODO: Get better estimate, and according to page 132, is effected bbu electomagnetic energy localized in a cavity
    lattice_constant = calculate_lattice_constant_in_units(wave_length, APPROXIMATLEY_LINEAR_FREQUENCY_REGION[0])
    radius = solve_for_radius(lattice_constant, radius_to_lattice_constant_ratio)
    print("Wavelength ", wave_length)
    print("Lattice Constant", lattice_constant)
    print("Radius ", radius)
    print("Mean Lifetime", mean_life_time_estimation)
    print("Speed of Light MICROMETERS", SPEED_OF_LIGHT_MICRO_METERS_PER_SECOND)
    top = gdspy.Cell("top")
    tk.add(top, OBLIGITORY_WAVE_GUIDE)
    single_stairwell_rods(
            top, 
            lattice_constant, 
            radius, 
            mean_life_time_estimation
        )
    print("Done generating geometry")
    tk.build_mask(top, OBLIGITORY_WAVE_GUIDE_TEMPLATE, final_layer = 200, final_datatype = 0)
    print("Done building mask, writing file")
    gdspy.write_gds('stairwell_mk0_static_with_gridded_wave_guide_first_mask.gds', unit=1.0e-6, precision=1.0e-9)
    print("Done writing file.")
    gdspy.LayoutViewer()

if __name__ == "__main__": 
   main()








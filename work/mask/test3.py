from enum import Enum
import numpy as np
import gdspy
from picwriter import toolkit as tk
import picwriter.components as pc
import custom_components as cc
#from custom_components.RectangularGrid2DTemplate import RectangleGrid2DTemplate#RectangularGrid2DTemplate#  as rg2dt
#from custom_components import RectangularGrid2D# as rg2d

UNITS = 1e-6 # Micrometers
NANOMETERS_IN_METERS = 1e-9
SPEED_OF_LIGHT_METERS_PER_SECOND = 299792458
SPEED_OF_LIGHT_MICRO_METERS_PER_SECOND = SPEED_OF_LIGHT_METERS_PER_SECOND / UNITS
MICROMETERS_IN_MILLIMETER = 10000
CHIP_TOTAL_WIDTH = MICROMETERS_IN_MILLIMETER * 5 # .5 centimeters

OBLIGITORY_WAVE_GUIDE_TEMPLATE = pc.WaveguideTemplate()
OBLIGITORY_WAVE_GUIDE = pc.Waveguide([(0, 1), (1, 1)], OBLIGITORY_WAVE_GUIDE_TEMPLATE)

APPROXIMATLEY_LINEAR_FREQUENCY_REGION = (.35, .42) # Photonic Crystal's Modeling the Flow of Light, Chapter 10 Page 194
                                                   # Wikipedia "Group Velocity" https://en.wikipedia.org/wiki/Group_velocity
                                                   # Various papers doubting superluminal phenomonia/the Hartman Effect citing the change in the shape of the wave packet.

def calculate_lattice_constant(
            wavelength : float, 
            frequency : float, 
            speed_of_light : float = SPEED_OF_LIGHT_MICRO_METERS_PER_SECOND
        ) -> float:
    return ((frequency / UNITS) * 2 * np.pi * speed_of_light) / wavelength

def solve_for_radius(lattice_constant : float, radius_over_lattice_constant_ratio : float) -> float: 
    return lattice_constant * radius_over_lattice_constant_ratio

def single_stairwell_rods(
            top, 
            lattice_constant : float, 
            radius : float, 
            mean_life_time : float, 
            speed_of_light = SPEED_OF_LIGHT_MICRO_METERS_PER_SECOND, 
            potential_height_ratios = [1, 2 / 3, 1 / 3], 
            segment_length_ratios = [1 / 3, 1 / 3, 1 / 3], 
            total_length = CHIP_TOTAL_WIDTH, 
            width = 7
        ): 
    grid_template = cc.RectangleGrid2DTemplate(lattice_constant, radius, radius)
    segment_unit_width = (mean_life_time * speed_of_light) / lattice_constant
    print(CHIP_TOTAL_WIDTH, lattice_constant)
    length_indicies = CHIP_TOTAL_WIDTH / lattice_constant
    print(length_indicies, width )
    wave_guide_grid = cc.RectangularGrid2D(
            grid_template, 
            (int(length_indicies), width), 
            [(int(np.floor(width / 2) + 1), -1)]
        )
    tk.add(top, waveGuideGrid)

def main(): 
    wave_length = 1550 * (NANOMETERS_IN_METERS / UNITS)
    radius_to_lattice_constant_ratio = .2 # TODO: Update from .2 for concodinanite. Photonic Crystal's Modeling the Flow of Light, Chapter 10 Page 192
    mean_life_time_estimation = 10e-12 / UNITS # TODO: Get better estimate, and according to page 132, is effected bbu electomagnetic energy localized in a cavity
    lattice_constant = calculate_lattice_constant(wave_length, APPROXIMATLEY_LINEAR_FREQUENCY_REGION[0])
    radius = solve_for_radius(lattice_constant, radius_to_lattice_constant_ratio)
    top = gdspy.Cell("top")
    tk.add(top, OBLIGITORY_WAVE_GUIDE)
    single_stairwell_rods(
            top, 
            lattice_constant, 
            radius, 
            mean_life_time_estimation
        )
    tk.build_mask(top, OBLIGITORY_WAVE_GUIDE_TEMPLATE, final_layer = 1000, final_datatype = 0)
    gdspy.LayoutViewer()

if __name__ == "__main__": 
   main()








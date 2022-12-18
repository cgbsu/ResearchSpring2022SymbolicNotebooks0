import custom_components as cc
from stairwell_builder import *
from lifetime_scaled_stairwell import *

total_length = 10000
top = gdspy.Cell("top")
wave_guide_position = (0, 200)
pad_template = FixedStaggeredPadGroupTemplate(
        cc.StaggeredMetalTemplate(clad_width = 0), 
        0, 
        100, 
        [],
        10
    )
wave_length : float = 1.5
freuquency : float = .35
admission : float = .2 # radius divided by lattice_constant
mean_life_time : float = 10e-12
speed_of_light : float = 3e14
lifetime_to_bondpad_length_ratio : float = 3
waveguide_index_width_in_lattice : float = 7
lattice_wave_guide : bool = True
layers = cc.RectangleGrid2DTemplate.DEFAULT_LAYERS
grid_stairwell = scaled_to_lifetime_stairwell(
        total_length, 
        top, 
        wave_guide_position, 
        pad_template, 
        wave_length, 
        freuquency, 
        admission, 
        mean_life_time,  
        speed_of_light, 
        lifetime_to_bondpad_length_ratio, 
        waveguide_index_width_in_lattice, 
        lattice_wave_guide, 
        layers
    )
slit_stairwell = scaled_to_lifetime_stairwell(
        total_length, 
        top, 
        (0, 1000), 
        pad_template, 
        wave_length, 
        freuquency, 
        admission, 
        mean_life_time,  
        speed_of_light, 
        lifetime_to_bondpad_length_ratio, 
        waveguide_index_width_in_lattice, 
        False, 
        layers
    )
smaller_grid_stairwell = scaled_to_lifetime_stairwell(
        total_length, 
        top, 
        (0, 3000), 
        pad_template, 
        wave_length, 
        freuquency, 
        admission, 
        mean_life_time,  
        speed_of_light, 
        3 * wave_length * (1 / mean_life_time / speed_of_light), 
        waveguide_index_width_in_lattice, 
        lattice_wave_guide, 
        layers, 
        clad_width = 0
    )
gdspy.LayoutViewer()



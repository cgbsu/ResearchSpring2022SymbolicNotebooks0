import sys
from pathlib import Path
import gdspy
sys.path.insert(0, '..')
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from stairwell_builder import *
import custom_components as cc

def calculate_lattice_constant(wavelength_in_units : float, frequency : float) -> float:
    # wavelength * lattice_constant / (2 * pi * c) = f
    # --> f * 2 * pi * c / wavelength = lattice_constatnt
    return frequency * wavelength_in_units

def solve_for_radius(lattice_constant : float, radius_over_lattice_constant_ratio : float) -> float: 
    return lattice_constant * radius_over_lattice_constant_ratio

def scaled_to_lifetime_stairwell(
            total_length, 
            top, 
            wave_guide_position, 
            pad_template : FixedStaggeredPadGroupTemplate, 
            wave_length : float, 
            frequency : float, 
            admission : float, # radius divided by lattice_constant
            mean_life_time : float, 
            speed_of_light : float, 
            lifetime_to_bondpad_length_ratio : float, 
            waveguide_index_width_in_lattice : float, 
            lattice_wave_guide : bool = True, 
            layers = cc.RectangleGrid2DTemplate.DEFAULT_LAYERS, 
            clad_width = 10
        ): 
    lattice_constant : float = calculate_lattice_constant(
            wave_length, 
            frequency
        )
    radius : float = solve_for_radius(lattice_constant, admission)
    if lattice_wave_guide == True: 
        extent_y_index : int = waveguide_index_width_in_lattice
        extent_x_index = int(round(total_length / lattice_constant))
        defect_layer_index = int(np.floor(extent_y_index / 2))
        grid_template = cc.RectangleGrid2DTemplate(
                lattice_constant, 
                radius, 
                layers = layers
            )
        wave_guide = cc.RectangularGrid2D(
                grid_template, 
                (extent_x_index, extent_y_index), 
                [(-1, defect_layer_index)], 
                position = wave_guide_position 
            )
    else: 
        wave_guide_template = pc.WaveguideTemplate(
                wg_type = "slot", 
                wg_width = lattice_constant * waveguide_index_width_in_lattice, 
                clad_width = clad_width
            )
        wave_guide = pc.Waveguide([wave_guide_position, (total_length, wave_guide_position[1])], wave_guide_template)

    tk.add(top, wave_guide)
    return Stairwell(
            top, 
            wave_guide, 
            pad_template, 
            place_lifetime_scaled_static_pad, 
            total_length, 
            builder_arguments = {
                    "speed_of_light" : speed_of_light, 
                    "mean_life_time" : mean_life_time
                }
        )

if __name__ == "__main__": 
    total_length = 10000
    top = gdspy.Cell("top")
    wave_guide_position = (0, 200)
    pad_template = FixedStaggeredPadGroupTemplate(
            cc.StaggeredMetalTemplate(), 
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
    stairwell = scaled_to_lifetime_stairwell(
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
    second_stairwell = scaled_to_lifetime_stairwell(
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
    gdspy.LayoutViewer()



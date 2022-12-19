import sys
from pathlib import Path
sys.path.insert(0, '..')
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
import gdspy
from picwriter import toolkit as tk
import custom_components as cc
from stairwell_builder import *
from lifetime_scaled_stairwell import *

SPEED_OF_LIGHT_METERS_PER_SECOND = 299792458
MICROMETERS_IN_METERS = 1e6
MICROMETERS_IN_CENTIMETER = MICROMETERS_IN_METERS / 100

def place_static_one_well_static_test_series(
            top : gdspy.Cell, 
            metal_template : pc.MetalTemplate, 
            pad_template : FixedStaggeredPadGroupTemplate, 
            start_position : float = (0, 0), 
            total_length : float = MICROMETERS_IN_CENTIMETER, 
            wave_length : float = 1.5, 
            freuquency : float = .35, 
            admission : float = .2, # radius divided by lattice_constant
            mean_life_time : float = 10e-12, 
            speed_of_light : float \
                    = SPEED_OF_LIGHT_METERS_PER_SECOND * MICROMETERS_IN_METERS, 
            waveguide_index_width_in_lattice : float = 7, 
            lattice_wave_guides : bool = True, 
            slot_wave_guides : bool = True, 
            split_total_length : bool = True, 
            layers = cc.RectangleGrid2DTemplate.DEFAULT_LAYERS, 
            wave_length_scale_tests : bool = True, 
            bond_pad_dimensions = (
                    .1 * MICROMETERS_IN_CENTIMETER, 
                    .1 * MICROMETERS_IN_CENTIMETER
                ), 
            trace_length_ratio = 2, 
            steps = 3
        ): 
    total_stairwell_pad_length = total_stairwell_lifetime_scaled_pad_length(
            wave_length, 
            mean_life_time, 
            speed_of_light, 
            steps, 
            False
        )
    bond_pad_width = bond_pad_dimensions[cc.DimensionalIndex.X.value]
    bond_pad_height = bond_pad_dimensions[cc.DimensionalIndex.Y.value]
    position = start_position
    def make_top_pad():
        return pc.Bondpad(
                metal_template, 
                bond_pad_width, 
                bond_pad_height, 
                port = (
                        position[cc.DimensionalIndex.X.value], 
                        position[cc.DimensionalIndex.Y.value] \
                                + (bond_pad_height * trace_length_ratio) * 2 \
                                + pad_template.width * 2 + waveguide_index_width_in_lattice 
                    )
            )
    bottom_pad = pc.Bondpad(
            metal_template, 
            bond_pad_width, 
            bond_pad_height, 
            port = position
        )
    top_pad = make_top_pad()
    tk.add(top, top_pad)
    tk.add(top, bottom_pad)
    wave_length_scales = [
            (10 ** power) * wave_length \
            for power in range(
                    int(np.log10(np.floor(wave_length))), 
                    int(np.log10(np.floor(total_stairwell_pad_length)))
                )
        ] + [
                wave_length * steps, 
                ((2 / 3) * steps) * wave_length, 
                wave_length, 
                2 * wave_length / steps, 
                wave_length / (2 * steps), 
                wave_length / steps
            ]
    lifetime_scales = [steps, (2 / 3) * steps, steps / 3]
    calculate_wave_guide_position = lambda : ( \
            position[cc.DimensionalIndex.X.value], 
            position[cc.DimensionalIndex.Y.value] \
                    + (bond_pad_height * trace_length_ratio) \
                    + pad_template.width + (waveguide_index_width_in_lattice / 2) \
        )
    wave_guide_position = calculate_wave_guide_position()
    def place_with_scales(scales, wave_length_scaled, half): 
        nonlocal top_pad
        nonlocal bottom_pad
        nonlocal wave_guide_position
        nonlocal position
        [False] if wave_length_scale_tests else [False, True]
        is_lattice_wave_guide = []
        if lattice_wave_guides == True: 
            is_lattice_wave_guide.append(True)
        if slot_wave_guides == True: 
            is_lattice_wave_guide.append(False)
        for lattice_wave_guide in [False]:#is_lattice_wave_guide: 
            for scale in scales: 
                grid_stairwell = scaled_to_lifetime_stairwell(
                        total_length / 2 if half == True else total_length, 
                        top, 
                        wave_guide_position, 
                        pad_template, 
                        wave_length, 
                        freuquency, 
                        admission, 
                        mean_life_time,  
                        speed_of_light, 
                        scale, 
                        waveguide_index_width_in_lattice, 
                        [top_pad.portlist["output"]["port"]], 
                        [bottom_pad.portlist["output"]["port"]], 
                        metal_template, 
                        lattice_wave_guide, 
                        layers, 
                        scale_to_wave_length = wave_length_scaled
                    )
                position = top_pad.port
                wave_guide_position = calculate_wave_guide_position()
                bottom_pad = top_pad
                top_pad = make_top_pad()
                tk.add(top, top_pad)
    if wave_length_scale_tests == True: 
        place_with_scales(wave_length_scales, True, False)
    #print("Lifetime scales")
    #place_with_scales(lifetime_scales, False, False)

if __name__ == "__main__": 
    top = gdspy.Cell("top")
    metal_template = pc.MetalTemplate()
    pad_template = FixedStaggeredPadGroupTemplate(
            cc.StaggeredMetalTemplate(clad_width = 0), 
            0, 
            100, 
            [],
            10
        )
    place_static_one_well_static_test_series(top, metal_template, pad_template, lattice_wave_guides = False)
    gdspy.LayoutViewer()



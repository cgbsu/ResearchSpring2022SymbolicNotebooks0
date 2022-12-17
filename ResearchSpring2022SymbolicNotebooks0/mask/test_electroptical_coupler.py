import gdspy
from picwriter import toolkit as tk
import picwriter.components as pc
import picwriter.picsim as ps


# Define materials
epsSiO2 = 1.444**2
epsSi = 3.55**2
etch_stack = [(epsSiO2, 1.89), (epsSi, 0.07), (epsSiO2, 2.04)]
mstack = ps.MaterialStack(vsize=4.0, default_stack=etch_stack, name="Si waveguide")
waveguide_stack = [(epsSiO2, 1.89), (epsSi, 0.22), (epsSiO2, 1.89)]
clad_stack = [(epsSiO2, 1.89), (epsSi, 0.05), (epsSiO2, 2.06)]
mstack.addVStack(layer=1, datatype=0, stack=waveguide_stack)
mstack.addVStack(layer=2, datatype=0, stack=clad_stack)


# Define geometry

top = gdspy.Cell("top")

waveGuideTemplate = pc.WaveguideTemplate(
        wg_width = 0.45, 
        clad_width = 10, 
        bend_radius = 100, 
        resist = '+', 
        fab = "ETCH", 
        wg_layer = 1, 
        wg_datatype = 0, 
        clad_layer = 2, 
        clad_datatype = 0
    )

metalTemplate = pc.MetalTemplate()

paddedWaveGuide = pc.Waveguide([(0, 1000), (1000, 1000)], waveGuideTemplate)

inputPadRout = pc.MetalRoute([(500, 0), (500, 700)], metalTemplate)
outputPadRout = pc.MetalRoute([(500, 2000), (500, 1300)], metalTemplate)

inputPad = pc.Bondpad(metalTemplate, **inputPadRout.portlist["input"])
anodePad = pc.Bondpad(metalTemplate, **inputPadRout.portlist["output"])

outputPad = pc.Bondpad(metalTemplate, **outputPadRout.portlist["input"])
cathodePad = pc.Bondpad(metalTemplate, **outputPadRout.portlist["output"])

tk.add(top, paddedWaveGuide)
tk.add(top, inputPad)
tk.add(top, anodePad)
tk.add(top, inputPadRout)
tk.add(top, outputPad)
tk.add(top, cathodePad)
tk.add(top, outputPadRout)

tk.build_mask(top, waveGuideTemplate, final_layer = 3, final_datatype = 0)

gdspy.LayoutViewer()

# Simulate

ps.compute_mode(
        waveGuideTemplate, 
        mstack, 
        res = 128, 
        wavelength = 1.55, 
        sx = 3.0, 
        sy = 3.0, 
        plot_mode_number = 1, 
        polarization = "TE"
    )




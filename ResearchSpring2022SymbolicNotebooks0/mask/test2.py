import gdspy
from picwriter import toolkit as tk
import picwriter.components as pc

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



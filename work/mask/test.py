import gdspy
from picwriter import toolkit as tk
import picwriter.components as pc

top = gdspy.Cell("top")

top.add(gdspy.Rectangle((0, 0), (1000, 1000), layer = 100, datatype = 0))

top.add(gdspy.Round((500, 500), 200))#tolerence = 

braggWaveGuideTemplate = pc.WaveguideTemplate(
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

braggWaveGuide = pc.Waveguide([(25, 25), (975, 25), (975, 500), (25, 500), (25, 975), (975, 975)], braggWaveGuideTemplate)

braggReflector = pc.DBR(braggWaveGuideTemplate, 100, 10, 1, 0)

top.add(braggReflector)

tk.add(top, braggWaveGuideTemplate)

tk.build_mask(top, braggWaveGuideTemplate, final_layer = 3, final_datatype = 0)

gdspy.LayoutViewer()


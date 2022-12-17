import gdspy
from picwriter import toolkit as tk
import picwriter.components as pc
import custom_components as cc

def initialize(): 
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
    waveGuide = pc.Waveguide([(0, 1000), (1000, 1000)], waveGuideTemplate)
    return top, waveGuideTemplate, [waveGuide]



def buildMask(top, features, waveGuideTemplate): 
    for feature in features: 
        if hasattr(feature, "integrate"):
            feature.integrate(top)
        else: 
            tk.add(top, feature)
    tk.build_mask(top, waveGuideTemplate, final_layer = 3, final_datatype = 0)


def main(): 
    top, waveGuideTemplate, features = initialize()
    buildMask(top, features, waveGuideTemplate)
    gdspy.LayoutViewer()

if __name__ == "__main__": 
    main()



import gdspy
from picwriter import toolkit as tk
import picwriter.components as pc
import picwriter.picsim as ps
import sys

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

def createPortList(features): 
    ports = []
    for feature in features: 
        ports.append(feature.portlist["input"])
        ports.append(feature.portlist["output"])
    return ports

def simulate(simulationCell, features, ports, waveGuideTemplate): 
    # Build material stack
    epsSiO2 = 1.444**2
    epsSi = 3.55**2
    etch_stack = [(epsSiO2, 1.89), (epsSi, 0.07), (epsSiO2, 2.04)]
    mstack = ps.MaterialStack(vsize=4.0, default_stack=etch_stack, name="Si waveguide")
    waveguide_stack = [(epsSiO2, 1.89), (epsSi, 0.22), (epsSiO2, 1.89)]
    clad_stack = [(epsSiO2, 1.89), (epsSi, 0.05), (epsSiO2, 2.06)]
    mstack.addVStack(layer=1, datatype=0, stack=waveguide_stack)
    mstack.addVStack(layer=2, datatype=0, stack=clad_stack)
    
    ps.compute_transmission_spectra(
            simulationCell, 
            mstack, 
            waveGuideTemplate, 
            ports, 
            port_vcenter=0,
            port_height=1.5 * 0.22, 
            port_width=1.5 * waveGuideTemplate.wg_width, 
            dpml=0.5,
            res=20, 
            wl_center=1.55, 
            wl_span=0.6, 
            fields=True,
            norm=True, 
            parallel=False
        )

def main(): 
    top, waveGuideTemplate, features = initialize()
    buildMask(top, features, waveGuideTemplate)
    if sys.argv[1] == "view": 
        gdspy.LayoutViewer()
    else: 
        ports = createPortList(features)
        simulate(top, features, ports, waveGuideTemplate)

if __name__ == "__main__": 
    main()



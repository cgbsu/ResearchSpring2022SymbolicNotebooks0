from picwriter import toolkit as tk
# Based off of the MetalTemplate class from https://github.com/DerekK88/PICwriter/blob/master/picwriter/components/electrical.py

class StaggeredMetalTemplate:
    """Template for electrical wires that contains some standard information about the fabrication process and metal wire.
    Keyword Args:
       * **bend_radius** (float): Radius of curvature for bends in the metal route.  Defaults to zero.
       * **width** (float): Width of the metal route as shown on the mask.  Defaults to 20.
       * **clad_width** (float): Width of the cladding (region next to route, mainly used for positive-type photoresists + etching, or negative-type and liftoff).  Defaults to 20.
       * **resist** (string): Must be either '+' or '-'.  Specifies the type of photoresist used.  Defaults to `'+'`.
       * **fab** (string): If 'ETCH', then keeps resist as is, otherwise changes it from '+' to '-' (or vice versa).  This is mainly used to reverse the type of mask used if the fabrication type is 'LIFTOFF'.  Defaults to `'ETCH'`.
       * **metal_layer** (int): Layer type used for metal route.  Defaults to 11.
       * **metal_datatype** (int): Data type used for metal route.  Defaults to 0.
       * **clad_layer** (int): Layer type used for cladding.  Defaults to 12.
       * **clad_datatype** (int): Data type used for cladding.  Defaults to 0.
    """

    def __init__(
                self,
                bend_radius=0,
                width=20.0,
                clad_width=20.0,
                resist="+",
                fab="ETCH",
                metal_layers=[11, 12, 13],
                metal_datatype=0,
                clad_layers=[14, 15, 16],
                clad_datatype=0,
            ):
        assert len(metal_layers) == len(clad_layers)
        self.name = tk.getCellName(
            "MetalTemplate"
        )  # Each MetalTemplate is given a unique name
        self.width = width
        self.bend_radius = bend_radius
        self.clad_width = clad_width
        if resist != "+" and resist != "-":
            raise ValueError(
                "Warning, invalid input for type resist in " "MetalTemplate"
            )
        if fab == "ETCH":
            self.resist = resist  # default state assumes 'etching'
        else:  # reverse resist type if liftoff or something else
            self.resist = "+" if resist == "-" else "-"

        self.metal_layers = metal_layers
        self.metal_datatype = metal_datatype
        self.clad_layers = clad_layers
        self.clad_datatype = clad_datatype


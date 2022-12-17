from picwriter import toolkit as tk
import picwriter.components as pc

class BondPads: 
    METAL_TEMPLATE = pc.MetalTemplate()
    def __init__(self, from_, to, interediateRoute = [], template = METAL_TEMPLATE): 
        self.from_ = from_
        self.to = to
        self.interediateRoute = interediateRoute 
        self.template = template
        self.routePoints = [self.to] + self.interediateRoute + [self.from_]
        self.route = pc.MetalRoute(self.route, metalTemplate)
        self.inputPad = pc.Bondpad(self.template, **self.route.portlist["input"])
        self.outputPad = pc.Bondpad(self.template, **self.route.portlist["output"])
        self.portlist = self.route.portlist
    def integrate(self, top): 
        tk.add(top, self.routePad)
        tk.add(top, self.inputPad)
        tk.add(top, self.outputPad)


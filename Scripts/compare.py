# code for comparing the output of Regina and t3m


import t3m
from xml.sax import saxutils, make_parser
from xml.sax.handler import feature_namespaces


class FindSurfaceWeights(saxutils.DefaultHandler):
    def __init__(self):
        self.surface_lists = []
        self.in_surface_list = 0
        self.in_surface = 0
        
    def startElement(self, name, attrs):
        if name == 'packet' and attrs.get("type") == "Normal Surface List":
            self.in_surface_list = 1
            self.curr_surface_list = []
        
        if name == 'surface':
            self.in_surface = 1
            self.weight_string = ""
            self.num_weights = int(attrs.get("len", None))

    def characters(self, ch):
        if self.in_surface:
            self.weight_string  += ch

    def endElement(self, name):
        if name == 'packet' and self.in_surface_list:
            self.in_surface_list = 0
            self.curr_surface_list.sort()
            self.surface_lists.append(self.curr_surface_list)
            
        if name == 'surface':
            self.in_surface = 0
            L = [0,]*self.num_weights
            nums = map(int, self.weight_string.split())
            for i in range(0, len(nums), 2):
                L[nums[i]] = nums[i+1]

            self.curr_surface_list.append(L)


# The function below extracts all the normal surfaces
#  from the given uncompressed .rga file.  Note: If there
#  is any kind of nesting of NormalSurfaceLists, I think
# the parser will get confused.

def extract_normal_surface_lists(file):
    parser = make_parser()
    dh = FindSurfaceWeights()
    parser.setFeature(feature_namespaces, 0)
    parser.setContentHandler(dh)
    parser.parse(file)
    L = dh.surface_lists
    return L

def compute_surfaces(manifold_name):
    M = t3m.get_mcomplex(manifold_name)
    M.find_normal_surfaces()
    final_list =[]
    for S in M.NormalSurfaces:
        vec = [0]*(3 * len(M))
        for i in range(len(M)):
            pos = [2, 1, 0][S.Quadtypes[i]]
            vec[3*i  + pos] = S.Coefficients[i]
        final_list.append(vec)
    final_list.sort()
    return final_list



    


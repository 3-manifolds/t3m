import gzip, struct, os, types, re
from SnapPea import *
from mcomplex import *

link_directory = os.path.join('Tables','Links') 
table_directory = os.path.join('Tables','Knots')

#------------------- Support for SnapPea Censuses --------------------
# Takes a list where the ith element represents the gluing data
# for the ith tetrahedron:
#
#  ( [Neighbors], [Glueings] )
#
# and creates the corresponding Mcomplex.

def Mcomplex_from_data(gluing_data):
    num_tets = len(gluing_data)
    tets = map(lambda x: Tetrahedron(), range(num_tets))
    for i in range(num_tets):
        neighbors, perms = gluing_data[i]
        for k in range(4):
            tets[i].attach(TwoSubsimplices[k], tets[neighbors[k]], perms[k])
    return Mcomplex(tets)

# A manifold_list is a subscriptable object containing Mcomplexes
# representing orientable closed manifolds from the SnapPea Census.

class manifold_list:
    def __init__(self, census, fill=1):
        self.census = census
        self.fill = fill

    def __repr__(self):
        return self.census.__repr__()

    def __len__(self):
        return len(self.census)

    def __getitem__(self, i):
        if not (0 <= i < len(self.census)):
            raise IndexError, "Requested manifold is out of range."
        manifold = self.census[i]
        gluing_data = SnapPeaC.get_gluing_data(manifold.triangulation, self.fill)
        return Mcomplex_from_data(gluing_data)    

closed_orientable = manifold_list( ClosedCensus(), fill=1)
five_tet_cusped = manifold_list( CuspedCensus(), fill=0)
six_tet_cusped_orientable = manifold_list( CuspedCensus(6), fill=0) 
six_tet_cusped_nonorientable = manifold_list( CuspedCensus(6,'n'), fill=0) 
seven_tet_cusped_orientable = manifold_list( CuspedCensus(7), fill=0) 
seven_tet_cusped_nonorientable = manifold_list( CuspedCensus(7,'n'), fill=0) 

# ---------- Support for the Hoste-Thistlethwaite knot table ----------
#
# We store the Hoste-Thistlethwaite tables in two gzipped files:
# alternating.gz and nonalternating.gz.  An alternating knot is stored
# in an array of 8 bytes.  To recover the Dowker-Thistlethwaite
# encoding of an alternating knot with k crossings, double the first k
# nybbles of the array, after replacing 0's by 16's.  A non-alternating
# knot is stored in an array of 10 bytes.  For a k-crossing knot the
# first k nybbles describe the Dowker-Thistlethwaite encoding of the
# corresponding alternating knot. This is padded with a zero nybble if
# k is odd.  The next two bytes describe the signs for the DT encoding
# of the nonalternating knot, with the high order bit of the first
# byte representing the sign of the first crossing. For knots with
# less than 16 crossings each record is padded with 0's.  But these
# are adequately handled by gzip, with the result that the two gzipped
# files are somewhat smaller than the carefully packed files that are
# supplied with the knotscape program.


# These dictionaries are used in accessing the tables.  The key is the
# number of crossings, the value is the number of knots with that many
# crossings.

Alternating_numbers = { 3:1, 4:1, 5:2, 6:3, 7:7, 8:18, 9:41, 10:123, 11:367,
                        12:1288, 13:4878, 14:19536, 15:85263, 16:379799 }

Nonalternating_numbers = { 8:3, 9:8, 10:42, 11:185, 12:888, 13:5110,
                           14:27346, 15:168030, 16:1008906 }

Alternating_offsets = {}
offset = 0
for i in range(3,17):
    Alternating_offsets[i] = offset
    offset +=  Alternating_numbers[i]
Num_Alternating = offset

Nonalternating_offsets = {}
offset = 0
for i in range(8,17):
    Nonalternating_offsets[i] = offset
    offset += Nonalternating_numbers[i]
Num_Nonalternating = offset

# These are the gzipped files holding the tables.
Alternating_table = gzip.open(os.path.join(table_directory, 'alternating.gz') )
Nonalternating_table = gzip.open(os.path.join(table_directory, 'nonalternating.gz') )

def extract_HT_knot(record, crossings, alternation):
    DT=[]
    size = (1+crossings)/2
    for byte in record[:size]:
        first_nybble = (byte & 0xf0) >> 4
        if first_nybble == 0: first_nybble = 16
        DT.append(2*first_nybble)
        second_nybble = byte & 0x0f
        if second_nybble == 0: second_nybble = 16
        DT.append(2*second_nybble)
    if alternation == 'n':
        signs = record[-2]<<8 | record[-1]
        mask = 0x8000
        for i in range(crossings):
            if (signs & (mask >> i)) == 0:
                DT[i] = -DT[i]
    return DT[:crossings]

def get_HT_knot(crossings, alternation, index):
    size = (1 + crossings)/2
    index -= 1
    if ( alternation == 'a'
         and crossings in Alternating_numbers.keys()
         and 0 <= index < Alternating_numbers[crossings] ):
        offset = 8*(Alternating_offsets[crossings] +  index)
        Alternating_table.seek(offset)
        data = Alternating_table.read(size)
        record = struct.unpack('%dB'%size, data)
    elif ( alternation == 'n'
         and crossings in Nonalternating_numbers.keys()
         and 0 <= index < Nonalternating_numbers[crossings] ):
        offset = 10*(Nonalternating_offsets[crossings] +  index)
        Nonalternating_table.seek(offset)
        data = Nonalternating_table.read(size+2)
        record = struct.unpack('%dB'%(size+2), data)
    else:
        print """
        You have entered a Hoste-Thistlethwaite knot name with an
        inappropriate index or number of crossings."""
        raise ValueError
    return extract_HT_knot(record, crossings, alternation)

def get_HT_knot_by_index(alternation, index):
    DT=[]
    crossings = 16
    if alternation == 'a':
        for i in range(3,17):
            if Alternating_offsets[i] > index:
                crossings = i-1
                break
        Alternating_table.seek(8*index)
        size = (1 + crossings)/2
        data = Alternating_table.read(size)
        record = struct.unpack('%dB'%size, data)
    if alternation == 'n':
        for i in range(8,17):
            if Nonalternating_offsets[i] > index:
                crossings = i-1
                break
        Nonalternating_table.seek(10*index)
        size = (1 + crossings)/2
        data = Nonalternating_table.read(size+2)
        record = struct.unpack('%dB'%(size+2), data)
    return extract_HT_knot(record, crossings, alternation)

# An HT_knot_list is a subscriptable object containing Mcomplexes
# representing exteriors of knots in the Hoste-Thistlethwaite tables.

class HT_knot_list:
    def __init__(self, alternation):
        if alternation != 'a' and alternation != 'n':
            alternation = None
        self.alternation = alternation

    def __repr__(self):
        if self.alternation == 'a':
            return """
    Exteriors of alternating knots from the Hoste-Thistlethwaite table."""
        elif self.alternation == 'n':
            return """
    Exteriors of non-alternating knots from the Hoste-Thistlethwaite table."""
        else:
            return """
    Improperly initialized knot list."""

    def __len__(self):
        if self.alternation == 'a':
            return Num_Alternating
        elif self.alternation == 'n':
            return Num_Nonalternating
        else:
            return 0
        
    def __getitem__(self, i):
        if not (0 <= i < self.__len__() ):
            raise IndexError, "Requested knot exterior is out of range."
        DT = get_HT_knot_by_index(self.alternation, i)
        manifold = Triangulation(SnapPeaC.get_triangulation_from_DT(DT)) 
        gluing_data = SnapPeaC.get_gluing_data(manifold.triangulation, 0)
        return Mcomplex_from_data(gluing_data)    

alternating_knot_ext = HT_knot_list('a')
nonalternating_knot_ext = HT_knot_list('n')

#---------- Support for the kitchen sink! ----------
#
# The following function loads a manifold specified by name, according
# to the following conventions:  
#
#   1. Numbers in parens at the end mean do Dehn filling on the loaded
#   manifold, e.g. m125(1,2)(4,5) means do (1,2) filling on the first
#   cusp and (4,5) filling on the second cusp.
#
#   2. Names of the form m123, s123, v123, and so on refer to the
#   SnapPea Census manifolds.
#
#   3. Names of the form 4_1, 04_1, 4_01, 5^2_6, 6_4^7, etc, refer to
#   complements of links in Rolfsen's table.  Similary L20935.  l104001, etc.
#
#   4. Names of the form b++LLR, b+-llR, bo-RRL, bn+LRLR load the
#   correponding torus bundle.
#
#   5. Names of the form 11a17 or 12n345 refer to complements of knots in
#   the Hoste-Thisthlethwaite tables.
#
#  If one of the above rules does _not_ apply, it looks for a file
#  with the specified name in the current directory and looks in the
#  path given by the user variable SNAPPEA_MANIFOLD_DIRECTORY.

split_filling_info = re.compile("(.*?)((?:\([0-9 .+-]+,[0-9 .+-]+\))+)")
is_census_manifold = re.compile("([msvxy])([0-9]+)$")
is_torus_bundle = re.compile("b([+-no])([+-])([lLrR]+)$")
is_knot_complement = re.compile("(?P<crossings>[0-9]+)_(?P<index>[0-9]+)$")
is_link_complement1 = re.compile("(?P<crossings>[0-9]+)[\^](?P<components>[0-9]+)[_](?P<index>[0-9]+)$")
is_link_complement2 = re.compile("(?P<crossings>[0-9]+)[_](?P<index>[0-9]+)[\^](?P<components>[0-9]+)$")
is_link_complement3 = re.compile("[lL]([0-9]+)")
is_HT_knot = re.compile('(?P<crossings>[0-9]*)(?P<alternation>[an])(?P<index>[0-9]*)')

spec_dict = {'m': (5, 1),
             's': (6, 1),
             'v': (7, 1),
             'x': (6, 0),
             'y': (7, 0)}

def get_manifold(name):
    # get filling info, if any
    m = split_filling_info.match(name)
    if m:
        real_name = m.group(1)
        fillings = re.subn("\)\(", "),(", m.group(2))[0]
        fillings = eval( "[" + fillings + "]" )
    else:
        real_name = name
        fillings = ()

    triangulation = None

    # Check for a census manifold
    m = is_census_manifold.match(real_name)
    if m:
         spec, orientable = spec_dict[m.group(1)]
         triangulation = SnapPeaC.get_cusped_census_manifold (
             spec, orientable, int(m.group(2)) )
        
    if triangulation == None:
     # Check for a puntured torus bundle 
        m = is_torus_bundle.match(real_name)
        if m:
            LRstring = m.group(3).upper()
            negative_determinant = negative_trace = 0

            if m.group(1) == '-' or m.group(1) == 'n':
                negative_determinant = 1
            
            if m.group(2) == '+':
                negative_trace = 0
            else:
                negative_trace = 1

            triangulation =  SnapPeaC.easy_triangulate_punctured_torus_bundle(
                negative_determinant, negative_trace, LRstring);
    
    if triangulation == None:
    # Check for a link
        filename = None
        m = is_knot_complement.match(real_name)
        if m:
            filename = "L1%.2d%.3d" % (int(m.group("crossings")),
                                       int(m.group("index")))
        m = is_link_complement1.match(real_name)
        if m:
            filename = "L%.1d%.2d%.3d" % (int(m.group("components")),
                                          int(m.group("crossings")),
                                          int(m.group("index")))
        m = is_link_complement2.match(real_name)
        if m:
            filename = "L%.1d%.2d%.3d" % (int(m.group("components")),
                                          int(m.group("crossings")),
                                          int(m.group("index")))
        m = is_link_complement3.match(real_name)
        if m:
            filename = "L" + m.group(1)

        if filename:
            try:
                 triangulation = SnapPeaC.get_triangulation(
                     os.path.join(link_directory, filename) )
            except:
                raise IOError, "Requested link complement " + real_name + " not found."

    if triangulation == None:
    # Check for a Hoste-Thistlethwaite knot.
        m = is_HT_knot.match(real_name)
        if m:
            DT = get_HT_knot(int(m.group("crossings")),
                             m.group("alternation"),
                             int(m.group("index")))
            triangulation = (SnapPeaC.get_triangulation_from_DT(DT)) 

    if triangulation == None:
    # If all else fails, try to load a manifold from a file.
        try:
            locations = [os.curdir, os.environ["SNAPPEA_MANIFOLD_DIRECTORY"]]
        except KeyError:
            locations = [os.curdir]
        found = 0
        for location in locations:
            filename = os.path.join(location,real_name)
            if os.path.isfile(filename):
                triangulation = SnapPeaC.get_triangulation(filename)
                break;
            
    if triangulation == None:
    # Give up if we failed to locate a manifold.
        raise IOError, "Requested manifold %s not found." % real_name
        
    # Otherwise do the dehn filling.
    if len(fillings) > 0:
        num_cusps = SnapPeaC.get_num_cusps(triangulation) 
        if len(fillings) > num_cusps:
            raise ValueError, "More fillings requested than manifold has cusps"
        for i in range(len(fillings)):
            m, l = fillings[i]
            SnapPeaC.set_cusp_info(triangulation, i, m, l, 0)
        s = real_name
        for filling in fillings:
            s = s + "(%s,%s)" % filling
        SnapPeaC.set_triangulation_name(triangulation, s)
        gluing_data = SnapPeaC.get_gluing_data(triangulation, 1)

    else:
        gluing_data = SnapPeaC.get_gluing_data(triangulation, 0)

    return Mcomplex_from_data(gluing_data)

# List of what gets exported from this module.

__all__ = ('closed_orientable',
           'five_tet_cusped',
           'six_tet_cusped_orientable',
           'six_tet_cusped_nonorientable',
           'seven_tet_cusped_orientable',
           'seven_tet_cusped_nonorientable',
           'alternating_knot_ext',
           'nonalternating_knot_ext',
           'get_manifold')


import xml.etree.ElementTree as ET
from dataclasses import dataclass,field,make_dataclass
from subprocess import Popen, PIPE, STDOUT
from math import sin,cos,tan,atan,pi
import re
import os


'''
Description:

    This is a python wrapper for AVL

    So far the following features have been implemented:

        the classes:

        configuration
            header
            surfaces
                sections
                    control
                    design
        
        allow an input file for avl (filename.avl) to be generated, saved and imported from a vsp3 geometry

        the classes:
        
            case
                constraint
                sweep ---> this is a type of constraint that additionally has an entire list of values. the sweep class is used to run eg. multiple AoAs

        allow a run file for avl (filename.run) to be generated and saved

        
        the class:

            totalforce

        parses the avl ouput of the same name from a .txt file


        the run function runs the avl binary using a configuration object and multiple case objects as input
    
'''


#utils
#region
AVL_BINARY_PATH = 'avl/avl.exe'

# functions for reading output files using re
def readFloatAfter(expression: str, filestr: str):
    match = re.search(r'%s[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?' % (expression),filestr)
    return float(match.group(0).replace(expression,''))

def readStrAfter(expression: str, filestr: str):
    match = re.search(r'%s(.*)' % (expression),filestr)
    return float(match.group(0).replace(expression,''))


DEG2RAD = pi/180.0

#XML parser
def recparse(root:ET):
    if root:
        fields = dict()
        for child in root:
            contents = recparse(child)
            if child.tag in fields.keys():
                if isinstance(fields[child.tag], list):
                    fields[child.tag].append(contents)
                else:
                    value = fields[child.tag]
                    fields[child.tag] = [value,contents]
            else:
                fields[child.tag] = contents
        containerclass = make_dataclass(root.tag,fields.keys(),repr=True, eq=False)
        container = containerclass(*fields.values())
        return container
    else:
        if not len(root.attrib):
            return root.text
        else:
            if 'Value' in root.attrib.keys():
                return float(root.attrib['Value'])
            else:
                #should not occur
                return None


def rm(fname):
    if os.path.isfile(fname):
        os.remove(fname)
    return
#endregion

#AVL file classes
#region
#   lifting surface geometry
#   control surfaces
#   design parameters (twist)
#   airfoil properties

@dataclass
class header:
    #attr   dtype   default
    name:   str
    Sref:   float
    Cref:   float
    Bref:   float
    XYZref: tuple[float]
    Mach:   float   = 0.0
    iYsym:  int     = 0
    iZsym:  int     = 0
    Zsym:   float   = 0.0
    CDp:    float   = None

    def __str__(self):
        str_list = [
            '%s' % self.name,
            '#Mach',
            '%f' % self.Mach,
            '#IYsym\tiZsym\tZsym',
            '%d %d %f' % (self.iYsym, self.iZsym, self.Zsym),
            '#Sref\tCref\tBref',
            '%f %f %f' % (self.Sref, self.Cref, self.Bref),
            '#Xref\tYref\tZref',
            '%f %f %f' % self.XYZref
            ]
        if self.CDp is not None:
            str_list.extend(['#CDp','%f' % self.CDp])

        str_list.append('#')

        return '\n'.join(str_list)

    @classmethod
    def _fromvsp3Vehicle(cls,Vehicle):
        name = Vehicle.ParmContainer.Name

        Sref = 0.0
        Bref = 0.0
        Cref = 0.0

        XYZref = (0.0,0.0,0.0)


        for Geom in Vehicle.Geom:
            if Geom.GeomBase.TypeName == 'Wing':
                if Geom.ParmContainer.WingGeom.TotalArea > Sref:
                    Sref = Geom.ParmContainer.WingGeom.TotalArea
                    Bref = Geom.ParmContainer.WingGeom.TotalSpan
                    Cref = Geom.ParmContainer.WingGeom.TotalChord

        return cls(name,Sref,Bref,Cref,XYZref)

@dataclass
class design:
    #attr   dtype   default
    name:   str
    weight: float

    def __str__(self):
        str_list = [
            '#',
            'DESIGN',
            self.name + '%f' % self.weight
            ]
        return '\n'.join(str_list)

@dataclass
class control:
    #attr   dtype   default
    name:   str
    gain:   float
    Xhinge: float
    XYZhvec: tuple[float,float,float]
    SgnDup: float

    def __str__(self):
        str_list = [
            '#',
            'CONTROL',
            '#name\tgain\tXhinge\tXYZhvec\tSgnDup',
            '%s %f %f' %(self.name,self.gain,self.Xhinge) + ' %f %f %f' % self.XYZhvec + ' %f' %self.SgnDup
            ]
        return '\n'.join(str_list)

@dataclass
class section:
    #attr   dtype   default
    XYZle:  tuple[float,float,float]
    chord:  float
    ainc:   float   = 0.0
    Nspan:  int     = None
    Sspace: float   = None
    X1:     float   = None
    X2:     float   = None
    NACA:   str     = None
    Afile:  str     = None
    design: list[design] = None
    control:list[control] = None
    Claf:   float   = None

    def __str__(self):
        str_list = [
            '#'*50,
            'SECTION',
            '#Xle\tYle\tZle\tChord\tAinc\tNspan\tSspace',
            '#'
            ]

        if self.Nspan is not None and self.Sspace is not None:
            str_list.append('%f %f %f' % self.XYZle + ' %f %f' % (self.chord, self.ainc),' %d %f' % (self.Nspan, self.Sspace))
        else:
            str_list.append('%f %f %f' % self.XYZle + ' %f %f' % (self.chord, self.ainc))

        if self.NACA is not None:
            if self.X1 is not None and self.X2 is not None:
                str_list.extend([
                    '#',
                    'NACA\t %f %f' % (self.X1, self.X2),
                    self.NACA
                ])
            else:
                str_list.extend([
                    '#',
                    'NACA',
                    self.NACA
                ])
        elif self.Afile is not None:
            pass

        if self.design is not None:
            for dsgn in self.design:
                str_list.append(str(dsgn))

        if self.control is not None:
            for ctrl in self.control:
                str_list.append(str(ctrl))

        if self.Claf is not None:
            str_list.extend([
                '#',
                'CLAF',
                '%f' % self.Claf
            ])

        return '\n'.join(str_list)

    @classmethod
    def _fromvsp3XSec(cls,XSec,XYZle,ainc):
        chord = XSec.ParmContainer.XSec.Tip_Chord

        if XSec.XSec.XSecCurve.XSecCurve.Type == '7':
            #NACA PROFILE
            tovc = XSec.XSec.XSecCurve.ParmContainer.XSecCurve.ThickChord
            camber = XSec.XSec.XSecCurve.ParmContainer.XSecCurve.Camber
            camberloc = XSec.XSec.XSecCurve.ParmContainer.XSecCurve.CamberLoc
            nacastr = '%d%d%02d' % (round(camber*100),round(camberloc*10),round(tovc*100))
            return cls(XYZle,chord,ainc,NACA=nacastr)
        else:
            return cls(XYZle,chord,ainc)

@dataclass
class surface:
    #attr   dtype   default
    name:   str
    Nchord: int
    Cspace: float
    Nspan:  int
    Sspace: float
    sections: list[section]
    translate: tuple[float,float,float]
    component:int   = None
    Ydupl:  float   = None
    scale:  tuple[float,float,float]   = None
    dAinc:  float   = None
    no_wake: bool   = False
    no_able: bool   = False
    no_load: bool   = False

    def __str__(self):
        str_list = [
            '#'*50,
            'SURFACE',
            '%s' % self.name,
            '#Nchord\tCspace\tNspan\tSspace'
            ]

        if self.Nspan is not None and self.Sspace is not None:
            str_list.append('%d %f %d %f' % (self.Nchord, self.Cspace, self.Nspan, self.Sspace))
        else:
            str_list.append('%d %f' % (self.Nchord, self.Cspace))

        if self.component is not None:
            str_list.extend([
                'COMPONENT',
                '%d' % self.component,
                '#'
            ])

        if self.Ydupl is not None:
            str_list.extend([
                'YDUPLICATE',
                '%f' % self.Ydupl,
                '#'
            ])

        if self.scale is not None:
            str_list.extend([
                'SCALE',
                '%f %f %f' % self.scale,
                '#'
            ])

        if self.translate is not None:
            str_list.extend([
                'TRANSLATE',
                '%f %f %f' % self.translate,
                '#'
            ])

        if self.dAinc is not None:
            str_list.extend([
                'ANGLE',
                '%f' % self.dAinc,
                '#'
            ])

        if self.no_wake:
            str_list.extend(['NOWAKE','#'])
        if self.no_able:
            str_list.extend(['NOABLE','#'])
        if self.no_load:
            str_list.extend(['NOLOAD','#'])



        for section in self.sections:
            str_list.append(str(section))

        return '\n'.join(str_list)

    @classmethod
    def _fromvsp3Geom(cls,Geom,Ydupl):
        name = Geom.ParmContainer.Name

        Nchord = 9
        Cspace = 1.0
        Nspan = 22
        Sspace = 2.0

        sections = list()

        XFormParams = Geom.ParmContainer.XForm
        translate = [XFormParams.X_Location,XFormParams.Y_Location,XFormParams.Z_Location]

        Xrot = XFormParams.X_Rotation
        cx = cos(Xrot*DEG2RAD)
        sx = sin(Xrot*DEG2RAD)
        Yrot = XFormParams.Y_Rotation
        cy = cos(Yrot*DEG2RAD)
        sy = sin(Yrot*DEG2RAD)

        reldihedralflag = bool(int(Geom.ParmContainer.WingGeom.RelativeDihedralFlag))
        reltwistflag = bool(int(Geom.ParmContainer.WingGeom.RelativeTwistFlag))

        for idx,XSec in enumerate(Geom.WingGeom.XSecSurf.XSec):
            XSecParams = XSec.ParmContainer.XSec
            if idx == 0:
                Xle = 0.0
                Yle = 0.0
                Zle = 0.0
                ainc = XSecParams.Twist

                #correct shift of root le due to inclination of wing
                dx = (1-cos(ainc*DEG2RAD))*XSecParams.Twist_Location*XSecParams.Tip_Chord
                dy = 0
                dz = (sin(ainc*DEG2RAD))*XSecParams.Twist_Location*XSecParams.Tip_Chord

                translate[0] +=   cy*dx       + 0*dy      + sy*dz
                translate[1] +=   sx*sy*dx    + cx*dy     - sx*cy*dz
                translate[2] +=  -cx*sy*dx    + sx*dy     + cx*cy*dz

            else:
                #leading edge sweep
                if bool(int(XSecParams.Sweep_Location)):
                    sweeple = XSecParams.Sweep_Location
                else:
                    sweeple = atan(tan(XSecParams.Sweep*DEG2RAD) + ( 2.0*XSecParams.Sweep_Location)*(1-XSecParams.Taper)/(XSecParams.Aspect * ( 1 + XSecParams.Taper)))

                #absolute dihedral
                if reldihedralflag:
                    if idx == 1:
                        absoluteDihedral = XSecParams.Dihedral
                        previousDihedral = XSecParams.Dihedral
                    else:
                        absoluteDihedral = XSecParams.Dihedral + previousDihedral
                        previousDihedral = XSecParams.Dihedral
                else:
                    absoluteDihedral = XSecParams.Dihedral

                #relative twist
                if reltwistflag:
                    relativeTwist = XSecParams.Twist
                else:
                    if idx == 1:
                        relativeTwist = XSecParams.Twist
                        previousTwist = XSecParams.Twist
                    else:
                        relativeTwist = XSecParams.Twist-previousTwist
                        previousTwist = XSecParams.Twist

                dx = tan(sweeple)*XSecParams.Span
                dy = XSecParams.Span*cos(absoluteDihedral*DEG2RAD)
                dz = XSecParams.Span*sin(absoluteDihedral*DEG2RAD)

                Xle +=  cy*dx       + 0*dy      + sy*dz
                Yle +=  sx*sy*dx    + cx*dy     - sx*cy*dz
                Zle += -cx*sy*dx    + sx*dy     + cx*cy*dz

                ainc += relativeTwist

            sections.append(section._fromvsp3XSec(XSec,(Xle,Yle,Zle),ainc))

        return cls(name,Nchord,Cspace,Nspan,Sspace,sections,tuple(translate),Ydupl=Ydupl)

@dataclass
class configuration:
    #attr   dtype   default
    header: header
    surfaces: list[surface]

    def __str__(self):
        str_list = [str(self.header)]

        for surface in self.surfaces:
            str_list.append(str(surface))

        return '\n'.join(str_list)

    def tofile(self,fname):
        with open(fname,'w') as file:
            file.write(str(self))
        return

    @classmethod
    def fromfile(cls,fname):
        if os.path.splitext(fname)[1] == '.vsp3':
            return cls._fromvsp3file(fname)

    @classmethod
    def _fromvsp3file(cls,fname):
        xmlroot = ET.parse(fname).getroot()
        Vehicle = recparse(xmlroot.find("./Vehicle"))

        avlheader = header._fromvsp3Vehicle(Vehicle)

        surfaces = list()

        for Geom in Vehicle.Geom:
            type = Geom.GeomBase.TypeName



            #symmetry
            planarsymmetryflag = int(Geom.ParmContainer.Sym.Sym_Planar_Flag)
            axialsymmetryflag = int(Geom.ParmContainer.Sym.Sym_Axial_Flag)
            ancestorflag = bool(int(Geom.ParmContainer.Sym.Sym_Ancestor_Origin_Flag))
            ancestorindex = int(Geom.ParmContainer.Sym.Sym_Ancestor)
            attachflag = int(Geom.ParmContainer.Attach.Trans_Attach_Flag)

            if planarsymmetryflag == 0 and axialsymmetryflag == 0:
                Ydupl =  None
            elif planarsymmetryflag == 2:
                #XZ symmetry
                if ancestorflag:
                    #geometry is symmetric w.r.t the attachment point of the geometry indicated by attachflag
                    ancestorID = Geom.GeomBase.ParentID
                    if attachflag==0 or ancestorID == None or ancestorindex == 0:
                        #geometry is symmetric w.r.t the global origin
                        Ydupl = 0.0
                    elif attachflag==1:
                        #geometry is symmetric w.r.t the parent part
                        for others in Vehicle.Geom:
                            if others.ParmContainer.ID == ancestorID:
                                Ydupl = others.ParmContainer.XForm.Y_Location
                                break
                    else:
                        raise NotImplementedError
                else:
                    if ancestorindex == 0:
                        #geometry is symmetric w.r.t the global origin
                        Ydupl = 0.0
                    else:
                        #geometry is symmetric w.r.t the (ancestorindex -1)th parent
                        currentGeom = Geom
                        while ancestorindex>1:
                            for others in Vehicle.Geom:
                                if others.ParmContainer.ID == ancestorID:
                                    currentGeom = others
                                    ancestorID = others.GeomBase.ParentID
                                    ancestorindex -= 1
                                    break
                        Ydupl = currentGeom.ParmContainer.XForm.Y_Location
            else:
                raise NotImplementedError


            if type == 'Wing':

                surfaces.append(surface._fromvsp3Geom(Geom,Ydupl))


        return cls(avlheader,surfaces)

#endregion

#run file classes
#region
#   AoA,AoS
#   flight conditions (Mach, altitude)
#   trim constraints (alpha such that CL = target, etc)
#   weight & balance
@dataclass
class constraint:
    #attr       dtype   default  unit
    tag:        str
    name:       str
    value:      float

    def __str__(self):
        return '{:<12s}=  {:f}'.format(self.name,self.value)

@dataclass(init=False)
class sweep(constraint):
    #attr       dtype   default  unit
    values:      list[float]

    def __init__(self,tag:str,name:str,values:list):
        self.tag = tag
        self.name = name
        self.value = values[0]
        self.values = values
    
    def __len__(self):
        return len(self.values)

    def __getitem__(self,idx):
        if isinstance(self.values[idx],list):
            newsweep = sweep(self.tag,self.name,self.values[idx])
            return newsweep
        else:
            self.value = self.values[idx]
            return self
    @classmethod
    def fromconstraint(cls,constr:constraint):
        return cls(constr.tag,constr.name,[constr.value])

@dataclass
class case:
    #attr       dtype   default  unit
    number:     int
    name:       str

    alpha:      constraint = field(init=False)
    beta:       constraint = field(init=False)
    pbov2V:     constraint = field(init=False)
    qcov2V:     constraint = field(init=False)
    rbov2V:     constraint = field(init=False)

    D:          dict[constraint] = field(default_factory=dict)

    CL:         float   = 0.0
    CD0:        float   = 0.0
    bank:       float   = 0.0    #deg
    elevation:  float   = 0.0    #deg
    heading:    float   = 0.0    #deg
    Mach:       float   = 0.0
    velocity:   float   = 0.0    #Lunit/Tunit
    density:    float   = 1.225  #Munit/Lunit^3
    grav_acc:   float   = 9.8067 #Lunit/Tunit^2
    turn_rad:   float   = 0.0    #Lunit
    load_fac:   float   = 1.0
    Xcg:        float   = 0.0    #Lunit
    Ycg:        float   = 0.0    #Lunit
    Zcg:        float   = 0.0    #Lunit
    mass:       float   = 1.0    #Munit
    Ixx:        float   = 1.0    #Munit-Lunit^2
    Iyy:        float   = 1.0    #Munit-Lunit^2
    Izz:        float   = 1.0    #Munit-Lunit^2
    Ixy:        float   = 1.0    #Munit-Lunit^2
    Iyz:        float   = 1.0    #Munit-Lunit^2
    Izx:        float   = 1.0    #Munit-Lunit^2
    dCLa:       float   = 0.0
    dCLu:       float   = 0.0
    dCMa:       float   = 0.0
    dCMu:       float   = 0.0

    def __post_init__(self):
        self.alpha = constraint('A','alpha',0.0)
        self.beta  = constraint('B','beta',0.0)
        self.pbov2V= constraint('R','pb/2V',0.0)
        self.qcov2V= constraint('P','qc/2V',0.0)
        self.rbov2V= constraint('Y','rb/2V',0.0)
        #self.CL    = constraint('C','CL',0.0)
        #self.CY    = constraint('S','CY',0.0)
        #self.Cl    = constraint('RM','Cl roll mom',0.0)
        #self.Cm    = constraint('PM','Cm pitchmom',0.0)
        #self.Cn    = constraint('YM','Cn yaw  mom',0.0)

    def __str__(self):
        str_list = [
            ' '+'-'*45,
            ' Run case' +'{:>3}:  {:}'.format(self.number,self.name),
            '',
            ' alpha        ->  ' + str(self.alpha),
            ' beta         ->  ' + str(self.beta),
            ' pb/2V        ->  ' + str(self.pbov2V),
            ' qc/2V        ->  ' + str(self.qcov2V),
            ' rb/2V        ->  ' + str(self.rbov2V)
            ]

        for key in self.D.keys():
            str_list.append(' {:<13}->  '.format(key) + str(self.D[key]))

        str_list += [
            ' alpha     =   %f     deg' % (self.alpha.value),
            ' beta      =   %f     deg' % (self.beta.value),
            ' pb/2V     =   %f'         % (self.pbov2V.value),
            ' qc/2V     =   %f'         % (self.qcov2V.value),
            ' rb/2V     =   %f'         % (self.rbov2V.value),
            ' CL        =   %f'         % (self.CL),
            ' CDo       =   %f'         % (self.CD0),
            ' bank      =   %f     deg' % (self.bank),
            ' elevation =   %f     deg' % (self.elevation),
            ' heading   =   %f     deg' % (self.heading),
            ' Mach      =   %f'         % (self.Mach),
            ' velocity  =   %f     m/s' % (self.velocity),
            ' density   =   %f     kg/m^3' % (self.density),
            ' grav.acc. =   %f     m/s^2'  % (self.grav_acc),
            ' turn_rad. =   %f     m'   % (self.turn_rad),
            ' load_fac. =   %f'         % (self.load_fac),
            ' X_cg      =   %f     m'   % (self.Xcg),
            ' Y_cg      =   %f     m'   % (self.Ycg),
            ' Z_cg      =   %f     m'   % (self.Zcg),
            ' mass      =   %f     kg'  % (self.mass),
            ' Ixx       =   %f     kgm^2' % (self.Ixx),
            ' Iyy       =   %f     kgm^2' % (self.Iyy),
            ' Izz       =   %f     kgm^2' % (self.Izz),
            ' Ixy       =   %f     kgm^2' % (self.Ixy),
            ' Iyz       =   %f     kgm^2' % (self.Iyz),
            ' Izx       =   %f     kgm^2' % (self.Izx),
            ' visc CL_a =   %f'         % (self.dCLa),
            ' visc CL_u =   %f'         % (self.dCLu),
            ' visc CM_a =   %f'         % (self.dCMa),
            ' visc CM_u =   %f'         % (self.dCMu),
            ' \n'
            ]

        return '\n'.join(str_list)

    def tofile(self,fname,append = False):
        if append:
            with open(fname,'a') as file:
                file.write(str(self))
            return
        else:
            with open(fname,'w') as file:
                file.write(str(self))
            return

    #def updatefromfile(self,fname):
    #
    #    with open(fname,'r') as file:
    #        lines = file.readlines()
    #        for i in range((len(lines)+1)//40):
    #            block = ''.join(lines[i*40:(i+1)*40])
    #            block = re.sub(' +',' ',block) #remove repeated whitespaces
    #            casenumber = readIntAfter('Run case ',block)
    #            casename = readStrAfter('Run case '+str(casenumber)+': ',block)
    #            if casenumber == self.number and casename == self.name:
    #
    #
    #                
    #
    #    return






#endregion

#output file classes
#region

@dataclass
class totalforces:
    conf: configuration
    case: case
    Alpha:   float = None
    Beta:    float = None
    Mach:    float = None
    pbov2V:  float = None
    qcov2V:  float = None
    rbov2V:  float = None
    pbov2Vst:float = None
    qcov2Vst:float = None
    rbpv2Vst:float = None
    CXtot:   float = None
    CYtot:   float = None
    CZtot:   float = None
    Cltot:   float = None
    Cmtot:   float = None
    Cntot:   float = None
    Cltotst: float = None
    Cmtotst: float = None
    Cntotst: float = None
    CLtot:   float = None
    CDtot:   float = None
    CDvis:   float = None
    CDind:   float = None
    CLff:    float = None
    CDff:    float = None
    CYff:    float = None
    e:       float = None

    @classmethod
    def fromfile(cls,fname,configuration,case):
        try:
            with open(fname,'r') as file:
                tf = cls(configuration,case)


                filestr = file.read()
                filestr = re.sub(' +',' ',filestr) #remove repeated whitespaces
                tf.Alpha = readFloatAfter('Alpha = ',filestr)
                tf.Beta  = readFloatAfter('Beta = ',filestr)
                tf.Mach  = readFloatAfter('Mach = ',filestr)
                tf.pbov2V= readFloatAfter('pb/2V = ',filestr)
                tf.qcov2V= readFloatAfter('qc/2V = ',filestr)
                tf.rbov2V= readFloatAfter('rb/2V = ',filestr)
                tf.pbov2Vst = readFloatAfter("p'b/2V = ",filestr)
                tf.qcov2Vst = tf.qcov2V
                tf.rbpv2Vst = readFloatAfter("r'b/2V = ",filestr)

                tf.CXtot = readFloatAfter('CXtot = ',filestr)
                tf.CYtot = readFloatAfter('CYtot = ',filestr)
                tf.CZtot = readFloatAfter('CZtot = ',filestr)
                tf.Cltot = readFloatAfter('Cltot = ',filestr)
                tf.Cmtot = readFloatAfter('Cmtot = ',filestr)
                tf.Cntot = readFloatAfter('Cntot = ',filestr)
                tf.Cltotst = readFloatAfter("Cl'tot = ",filestr)
                tf.Cmtotst = tf.Cmtot
                tf.Cntotst = readFloatAfter("Cn'tot = ",filestr)
                tf.CLtot = readFloatAfter('CLtot = ',filestr)
                tf.CDtot = readFloatAfter('CDtot = ',filestr)
                tf.CDvis = readFloatAfter('CDvis = ',filestr)
                tf.CDind = readFloatAfter('CDind = ',filestr)
                tf.CLff  = readFloatAfter('CLff = ',filestr)
                tf.CDff  = readFloatAfter('CDff = ',filestr)
                tf.CYff  = readFloatAfter('CYff = ',filestr)
                tf.e     = readFloatAfter('e = ',filestr)

            return tf
        except FileNotFoundError:
            return False



#endregion

# driver code
#region
#   execute the AVL binary
def run(config: configuration, cases: list[case], binarypath: str = AVL_BINARY_PATH):
    proc = Popen([binarypath],stdout=PIPE,stdin=PIPE,stderr=STDOUT,encoding='utf-8')

    # write configuration to .avl file
    conffilename = config.header.name + '.avl'
    rm(conffilename)
    config.tofile(conffilename)
    
    #write cases to .run file
    casefilename = config.header.name + '.run'
    rm(casefilename)
    cases[0].tofile(casefilename)
    for cs in cases[1:]:
        cs.tofile(casefilename,append=True)



    # define inputs for avl
    cmd = ""
    cmd += 'plop\nG\n\n'                        # turn off plotting
    cmd += 'load %s\n' % (conffilename)         # load configuration
    cmd += 'oper\n'                             # enter OPER mode
    '''------------------OPER------------------'''
    
    cmd += 'f\n%s\n' % (casefilename)           # load casefile
    for cs in cases:
        cmd += '{}\n'.format(cs.number)         #select case

        if not isinstance(cs.alpha,sweep):
            cs.alpha = sweep.fromconstraint(cs.alpha)
        if not isinstance(cs.beta,sweep):
            cs.beta = sweep.fromconstraint(cs.beta)
        if not isinstance(cs.pbov2V,sweep):
            cs.pbov2V = sweep.fromconstraint(cs.pbov2V)
        if not isinstance(cs.qcov2V,sweep):
            cs.qcov2V = sweep.fromconstraint(cs.qcov2V)
        if not isinstance(cs.rbov2V,sweep):
            cs.rbov2V = sweep.fromconstraint(cs.rbov2V)

        for i in range(len(cs.alpha)):
            for j in range(len(cs.beta)):
                for k in range(len(cs.pbov2V)):
                    for l in range(len(cs.qcov2V)):
                        for m in range(len(cs.rbov2V)):
                            tfname = 'tf_%s%d_%s%d_%s%d_%s%d_%s%d.txt' % (cs.alpha.tag,i,cs.beta.tag,j,cs.pbov2V.tag,k,cs.qcov2V.tag,l,cs.rbov2V.tag,m)
                            rm(tfname)
                            cmd += 'a %s %f\n' % (cs.alpha.tag,cs.alpha.values[i])      #set constraint for alpha
                            cmd += 'b %s %f\n' % (cs.beta.tag,cs.beta.values[j])        #set constraint for beta
                            cmd += 'r %s %f\n' % (cs.pbov2V.tag,cs.pbov2V.values[k])    #set constraint for roll rate
                            cmd += 'p %s %f\n' % (cs.qcov2V.tag,cs.qcov2V.values[l])    #set constraint for pitch rate
                            cmd += 'y %s %f\n' % (cs.rbov2V.tag,cs.rbov2V.values[m])    #set constraint for yaw rate

                            cmd += 'x\n'                    #run case
                            cmd += 'ft\n%s\n' % (tfname)    #save total forces 


    cmd += '\n'                                 #exit from oper menu
    '''----------------------------------------'''
    cmd += 'q\n'                                #quit

    output = proc.communicate(cmd)[0]

    #print(output)

    tfout = list()

    cs = cases[0]
    for i in range(len(cs.alpha)):
        for j in range(len(cs.beta)):
            for k in range(len(cs.pbov2V)):
                for l in range(len(cs.qcov2V)):
                    for m in range(len(cs.rbov2V)):
                        tfname = 'tf_%s%d_%s%d_%s%d_%s%d_%s%d.txt' % (cs.alpha.tag,i,cs.beta.tag,j,cs.pbov2V.tag,k,cs.qcov2V.tag,l,cs.rbov2V.tag,m)
                        forces = totalforces.fromfile(tfname,config,cs)
                        tfout.append(forces)
                        rm(tfname)

    rm(conffilename)
    rm(casefilename)

    return tfout



#todo
# sweep one single constraint (eg. sweep alpha -> CL = 0.7... 0.9   use result for CL/CD plot)



#endregion


if __name__ == '__main__':

    # example case: define avl configuration consisting of a single wing with non-zero incidence and twist
    #               define a case that trims the "aircraft" in pitch & run it

    avlheader = header('testcase',Sref=10,Cref=1,Bref=10,XYZref=(0.0,0.0,0.0))
    rootsec = section((0.0,0.0,0.0),1.0,ainc = 2.0)
    tipsec = section((0.5,5.0,0.0),0.5,ainc= 4.0)
    mainwing = surface('main wing',8,1.0,20,2.0,[rootsec,tipsec],(0.0,0.0,0.0),Ydupl = 0.0)
    testconfig = configuration(avlheader,[mainwing])
    #testconfig.tofile('testcase.avl')

    testcase = case(1,'test')
    #testcase.tofile('testcase.run')

    testcase.alpha = constraint('PM','Cm pitchmom',0.0)

    tfo = run(testconfig,[testcase])[0]

    print('configuration: %s' %tfo.conf.header.name)
    print('AoA          : %f' %tfo.Alpha)
    print('CL           : %f' %tfo.CLtot)
    print('Pitch mom    : %f' %tfo.Cmtot)




from dataclasses import dataclass,field
from enum import Enum
from re import compile,sub,split,findall
import math

from matplotlib.pyplot import close


@dataclass
class header:
    #attr   dtype   default
    name:   str
    Sref:   float
    Cref:   float
    Bref:   float
    XYZref: tuple[float,float,float]
    Mach:   float   = 0.0
    iYsym:  int     = 0
    iZsym:  int     = 0
    Zsym:   float   = 0.0
    CDp:    float   = None

    @property
    def configfilestring(self) -> str:
        str_list = [
            '%s' % self.name,
            '#Mach',
            '%f' % self.Mach,
            '#iYsym\tiZsym\tZsym',
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
    def _fromstring(cls,string):
        lines = string.strip().split('\n')
        name = lines[0]
        temp = lines[1].split()
        Mach = float(temp[0])
        temp = lines[2].split()
        iYsym,iZsym,Zsym = int(float(temp[0])),int(float(temp[1])),float(temp[2])
        temp = lines[3].split()
        Sref,Cref,Bref = float(temp[0]),float(temp[1]),float(temp[2])
        temp = lines[4].split()
        XYZref = (float(temp[0]),float(temp[1]),float(temp[2]))
        if len(lines)>5:
            temp = lines[5].split()
            CDp = float(temp[0])
        else:
            CDp = None
        return cls(name,Sref,Cref,Bref,XYZref,Mach,iYsym,iZsym,Zsym,CDp)

@dataclass
class design:
    #attr   dtype   default
    name:   str
    weight: float = 0.0

    @property
    def configfilestring(self) -> str:
        str_list = [
            '#',
            'DESIGN',
            self.name + '%f' % self.weight
            ]
        return '\n'.join(str_list)

    @classmethod
    def _fromstring(cls,string):
        temp = string.split()
        if len(temp) == 2:
            name,weight = temp[0],float(temp[1])
            return cls(name,weight)
        else:
            name = temp[0]
            return cls(name)
    
    @classmethod
    def lerp(cls,u,a,b):
        _lerp = lambda u,a,b: (1-u)*a + u*b
        name = a.name
        weight = _lerp(u,a.weight,b.weight)
        return cls(name,weight)
        
@dataclass
class control:
    #attr   dtype   default
    name:   str
    gain:   float
    Xhinge: float
    XYZhvec: tuple[float,float,float]
    SgnDup: float = None

    @property
    def configfilestring(self) -> str:
        str_list = [
            '#',
            'CONTROL']
        if self.SgnDup is not None:
            str_list.extend(['#name\tgain\tXhinge\tXYZhvec\tSgnDup',
            '%s %f %f' %(self.name,self.gain,self.Xhinge) + ' %f %f %f' % self.XYZhvec + ' %f' %self.SgnDup])
        else:
            str_list.extend(['#name\tgain\tXhinge\tXYZhvec\tSgnDup',
            '%s %f %f' %(self.name,self.gain,self.Xhinge) + ' %f %f %f' % self.XYZhvec])
        return '\n'.join(str_list)
    
    @classmethod
    def _fromstring(cls,string):
        temp = string.split()
        if len(temp)==7:
            name,gain,Xhinge,XYZhvec,SgnDup = temp[0],float(temp[1]),float(temp[2]),(float(temp[3]),float(temp[4]),float(temp[5])),float(temp[6])
            return cls(name,gain,Xhinge,XYZhvec,SgnDup)
        elif len(temp)==6:
            name,gain,Xhinge,XYZhvec = temp[0],float(temp[1]),float(temp[2]),(float(temp[3]),float(temp[4]),float(temp[5]))
            return cls(name,gain,Xhinge,XYZhvec)


    @classmethod
    def lerp(cls,u,a,b):
        _lerp = lambda u,a,b: (1-u)*a + u*b
        name = a.name
        gain = _lerp(u,a.gain,b.gain)
        Xhinge = _lerp(u,a.Xhinge,b.Xhinge)

        XYZhvec = ( _lerp(u,a.XYZhvec[0],b.XYZhvec[0]),
                    _lerp(u,a.XYZhvec[1],b.XYZhvec[1]),
                    _lerp(u,a.XYZhvec[2],b.XYZhvec[2]))
        
        SgnDup = a.SgnDup
        return cls(name,gain,Xhinge,XYZhvec,SgnDup)

class aftypes(Enum):
    NACA = 0
    AFILE = 1
    AIRFOIL = 2

@dataclass
class airfoil:
    type: aftypes
    X1:     float
    X2:     float
    NACA: str
    Afile: str
    Airfoil: tuple[tuple[float,float]]
    Claf: float
    CDCL: tuple[float,float,float,float,float,float]

    @property
    def configfilestring(self)->str:
        if self.type is aftypes.NACA:
            if self.X1 is not None and self.X2 is not None:
                str_list = [
                    '#',
                    'NACA\t %f %f' % (self.X1, self.X2),
                    self.NACA]
            else:
                str_list = [
                    '#',
                    'NACA',
                    self.NACA]
            return '\n'.join(str_list)
        elif self.type is aftypes.AFILE:
            if self.X1 is not None and self.X2 is not None:
                str_list = [
                    '#',
                    'AFILE\t %f %f' % (self.X1, self.X2),
                    self.Afile]
            else:
                str_list = [
                    '#',
                    'AFILE',
                    self.Afile]
            return '\n'.join(str_list)
        elif self.type is aftypes.AIRFOIL:
            if self.X1 is not None and self.X2 is not None:
                str_list = [
                    '#',
                    'AIRFOIL\t %f %f' % (self.X1, self.X2)]
                for x,y in self.Airfoil:
                    str_list.append('%f %f' % (x,y))
            else:
                str_list = [
                    '#',
                    'AIRFOIL']
                for x,y in self.Airfoil:
                    str_list.append('%f %f' % (x,y))
                str_list.append('%f %f' % (x,y))
            return '\n'.join(str_list)
    
    @classmethod
    def _fromstring(cls,string):
        lines = string.split('\n')
        type = None
        NACA,Afile,Airfoil,X1,X2 = None,None,None,None,None
        Claf = None
        CDCL = None

        while lines:
            line = lines.pop(0)
            keyw = line[0:4]
            if keyw == 'NACA':
                type = aftypes.NACA
                line = lines.pop(0)
                temp = line.split()
                if len(temp) >= 3:
                    NACA,X1,X2 = temp[0],float(temp[1]),float(temp[2])
                else:
                    NACA = temp[0]
            elif keyw == 'AFIL':
                type = aftypes.AFILE
                line = lines.pop(0)
                temp = line.split()
                if len(temp) >= 3:
                    Afile,X1,X2 = temp[0],float(temp[1]),float(temp[2])
                else:
                    Afile = temp[0]
            elif keyw == 'AIRF':
                type = aftypes.AIRFOIL
                temp = line.split()
                if len(temp) >= 3:
                    X1,X2 = float(temp[1]),float(temp[2])
                coords = list()
                while lines:
                    try:
                        temp = lines[0].split()
                        x,y = float(temp[0]),float(temp[1])
                        coords.append((x,y))
                        lines.pop(0)
                    except IndexError:
                        break
                Airfoil = tuple(coords)
            elif keyw == 'CLAF':
                line = lines.pop(0)
                temp = line.split()
                Claf = float(temp[0])
            elif keyw == 'CDCL':
                line = lines.pop(0)
                temp = line.split()
                CDCL = (float(temp[0]),float(temp[1]),float(temp[2]),float(temp[3]),float(temp[4]),float(temp[5]))

        if type == None:
            return None
        else:
            return cls(type,X1,X2,NACA,Afile,Airfoil,Claf,CDCL)

    @classmethod
    def lerp(cls,u,a,b):
        _lerp = lambda u,a,b: (1-u)*a + u*b
        
        type = None
        NACA,Afile,Airfoil,X1,X2 = None,None,None,None,None
        Claf = None
        CDCL = None

        if not (None in [a.X1,a.X2,b.X1,b.X2]):
            X1 = _lerp(u,a.X1,b.X1)
            X2 = _lerp(u,a.X2,b.X2)

        if not (None in [a.NACA,b.NACA]):
            camber      = round(_lerp(u,float(a.NACA[0]),float(b.NACA[0])))
            camberloc   = round(_lerp(u,float(a.NACA[1]),float(b.NACA[1])))
            tovc        = round(_lerp(u,float(a.NACA[2:]),float(b.NACA[2:])))
            NACA = "%d%d%02d" % (camber,camberloc,tovc)
            type = aftypes.NACA

        if not (None in [a.Claf,b.Claf]):
            Claf = _lerp(u,a.Claf,b.Claf)
        
        if not (None in [a.CDCL,b.CDCL]):
            CLmin = _lerp(u,a.CDCL[0],b.CDCL[0])
            CDmin = _lerp(u,a.CDCL[1],b.CDCL[1])
            CL0 = _lerp(u,a.CDCL[2],b.CDCL[2])
            CD0 = _lerp(u,a.CDCL[3],b.CDCL[3])
            CLmax = _lerp(u,a.CDCL[4],b.CDCL[4])
            CDmax = _lerp(u,a.CDCL[5],b.CDCL[5])
            CDCL = (CLmin,CDmin,CL0,CD0,CLmax,CDmax)
        
        return cls(type,X1,X2,NACA,Afile,Airfoil,Claf,CDCL)

@dataclass
class section:
    #attr   dtype   default
    XYZle:  tuple[float,float,float]
    chord:  float
    ainc:   float   = 0.0
    Nspan:  int     = None
    Sspace: float   = None
    Airfoil: airfoil = None
    designs: list[design] = None
    controls:list[control] = None


    @property
    def configfilestring(self) -> str:
        str_list = [
            '!'+'-'*49,
            'SECTION',
            '#Xle\tYle\tZle\tChord\tAinc\tNspan\tSspace']

        if self.Nspan is not None and self.Sspace is not None:
            str_list.append('%f %f %f' % self.XYZle + ' %f %f' % (self.chord, self.ainc) + ' %d %f' % (self.Nspan, self.Sspace))
        else:
            str_list.append('%f %f %f' % self.XYZle + ' %f %f' % (self.chord, self.ainc))

        if self.Airfoil is not None:
            str_list.append(self.Airfoil.configfilestring)

        if self.designs is not None:
            for dsgn in self.designs:
                str_list.append(dsgn.configfilestring)


        if self.controls is not None:
            for ctrl in self.controls:
                str_list.append(ctrl.configfilestring)


        return '\n'.join(str_list)

    @classmethod
    def _fromstring(cls,string):
        CONT = compile('CONT.*?\n|DESI.*?\n')
        segments = split(CONT,string)
        types = findall(CONT,string)

        Nspan,Sspace = None,None
        Airfoil = None
        designs = None
        controls = None

        
        temp = segments.pop(0).split('\n',1)
        airfoilstring = temp[1]
        temp = temp[0].split()
        if len(temp) == 5:
            XYZle,chord,ainc = (float(temp[0]),float(temp[1]),float(temp[2])),float(temp[3]),float(temp[4])
        elif len(temp) == 7:
            XYZle,chord,ainc,Nspan,Sspace = (float(temp[0]),float(temp[1]),float(temp[2])),float(temp[3]),float(temp[4]),int(float(temp[5])),float(temp[6])

        Airfoil = airfoil._fromstring(airfoilstring)
        controls = list()
        designs = list()
        for type,segment in zip(types,segments):
            if type[0:4] == 'CONT':
                controls.append(control._fromstring(segment))
            elif type[0:4] == 'DESI':
                designs.append(design._fromstring(segment))  
        return cls(XYZle,chord,ainc,Nspan,Sspace,Airfoil,designs,controls)

    @classmethod
    def lerp(cls,u,a,b):
        _lerp = lambda u,a,b: (1-u)*a + u*b
        XYZle = (_lerp(u,a.XYZle[0],b.XYZle[0]),
                 _lerp(u,a.XYZle[1],b.XYZle[1]),
                 _lerp(u,a.XYZle[2],b.XYZle[2]))
        chord = _lerp(u,a.chord,b.chord)
        ainc = _lerp(u,a.ainc,b.ainc)



        Nspan,Sspace = None,None
        Airfoil = None

        if not (None in [a.Nspan,a.Sspace,b.Nspan,b.Sspace]):
            Nspan = int(_lerp(u,float(a.Nspan),float(b.Nspan)))
            Sspace = _lerp(u,a.Sspace,b.Sspace)
        
        if not (None in [a.Airfoil,b.Airfoil]):
            Airfoil = airfoil.lerp(u,a.Airfoil,b.Airfoil)


        designs = None
        if not (None in [a.designs,b.designs]):
            designs = list()
            for desga in a.designs:
                for desgb in b.designs:
                    if desgb.name == desga.name:
                        designs.append(design.lerp(u,desga,desgb))

        controls = None
        if not (None in [a.controls,b.controls]):
            controls = list()
            for ctrla in a.controls:
                for ctrlb in b.controls:
                    if ctrlb.name == ctrla.name:
                        controls.append(control.lerp(u,ctrla,ctrlb))
    
           
        return cls(XYZle,chord,ainc,Nspan,Sspace,Airfoil,designs,controls)

@dataclass
class surface:
    #attr   dtype   default
    name:   str
    Nchord: int
    Cspace: float
    Nspan:  int
    Sspace: float
    sections: list[section]
    translate: tuple[float,float,float] = None
    component:int   = None
    Ydupl:  float   = None
    scale:  tuple[float,float,float]   = None
    dAinc:  float   = None
    no_wake: bool   = False
    no_albe: bool   = False
    no_load: bool   = False
    CDCL: tuple[float,float,float,float,float,float] = None

    @property
    def configfilestring(self) -> str:
        str_list = [
            '!'+'='*49,
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
        if self.no_albe:
            str_list.extend(['NOABLE','#'])
        if self.no_load:
            str_list.extend(['NOLOAD','#'])
        if self.CDCL is not None:
            str_list.extend([
                '#',
                'CDCL',
                '%f %f %f %f %f %f' % self.CDCL
            ])


        for section in self.sections:
            str_list.append(section.configfilestring)

        return '\n'.join(str_list)

    @classmethod
    def _fromstring(cls,string):
        SECT = compile('SECT.*?\n')
        segments = split(SECT,string)
        #read surface parameters
        lines = segments.pop(0).split('\n')
        name = lines.pop(0)
        temp = lines.pop(0).split()
        if len(temp) == 2:
            Nchord,Cspace,Nspan,Sspace = int(float(temp[0])),float(temp[1]),None,None
        elif len(temp) == 4:
            Nchord,Cspace,Nspan,Sspace = int(float(temp[0])),float(temp[1]),int(float(temp[2])),float(temp[3])
        component = None
        Ydupl = None
        scale = None
        translate = None
        dAinc = None
        no_wake = False
        no_albe = False
        no_load = False
        CDCL = None

        while lines:
            line = lines.pop(0)
            keyw = line[0:4]
            if keyw == 'COMP':
                line = lines.pop(0)
                temp = line.split()
                component = int(float(temp[0]))
            elif keyw == 'YDUP':
                line = lines.pop(0)
                temp = line.split()
                Ydupl = float(temp[0])
            elif keyw == 'SCAL':
                line = lines.pop(0)
                temp = line.split()
                scale = (float(temp[0]),float(temp[1]),float(temp[2]))
            elif keyw == 'TRAN':
                line = lines.pop(0)
                temp = line.split()
                translate = (float(temp[0]),float(temp[1]),float(temp[2]))
            elif keyw == 'ANGL':
                line = lines.pop(0)
                temp = line.split()
                dAinc = float(temp[0])
            elif keyw == 'NOWA':
                no_wake = True
            elif keyw == 'NOAL':
                no_albe = True
            elif keyw == 'NOLO':
                no_load =True
            elif keyw == 'CDCL':
                line = lines.pop(0)
                temp = line.split()
                CDCL = (float(temp[0]),float(temp[1]),float(temp[2]),float(temp[3]),float(temp[4]),float(temp[5]))
        
        sections = list()
        # read sections
        for segment in segments:
            sections.append(section._fromstring(segment))
        
        return cls(name,Nchord,Cspace,Nspan,Sspace,sections,translate,component,Ydupl,scale,dAinc,no_wake,no_albe,no_load,CDCL)         

    def splitspan(self,ys,zs):
        self.sections.sort(key= lambda x: math.sqrt(math.pow(x.XYZle[1],2) + math.pow(x.XYZle[2],2)))

        closestidx = min(range(len(self.sections)), key= lambda i: math.sqrt(math.pow(self.sections[i].XYZle[1] - ys,2) + math.pow(self.sections[i].XYZle[2] - zs,2)))

        radius_closest_section = math.sqrt(math.pow(self.sections[closestidx].XYZle[1],2) + math.pow(self.sections[closestidx].XYZle[2],2))

        radius_split_plane = math.sqrt(math.pow(ys,2) + math.pow(zs,2))

        if radius_closest_section > radius_split_plane:
            inboard_idx = closestidx - 1
            outboard_idx = closestidx
        else:
            inboard_idx = closestidx
            outboard_idx = closestidx + 1


        dYle = abs(self.sections[outboard_idx].XYZle[1] - self.sections[inboard_idx].XYZle[1])
        dZle = abs(self.sections[outboard_idx].XYZle[2] - self.sections[inboard_idx].XYZle[2])
        
        if dYle > 0:
            u = abs(ys - self.sections[outboard_idx].XYZle[1])/dYle
        elif dZle > 0:
            u = abs(zs - self.sections[outboard_idx].XYZle[2])/dZle
        else:
            raise Exception("Duplicated section")

        self.sections.insert(outboard_idx, section.lerp(u,self.sections[outboard_idx],self.sections[inboard_idx]))
        return


@dataclass
class body:
    #attr   dtype   default
    name:   str
    Nbody:  int
    Bspace: float
    translate: tuple[float,float,float]
    Bfile:  str
    Ydupl:  float   = None
    scale:  tuple[float,float,float]   = None
    X1:     float   = None
    X2:     float   = None

    @property
    def configfilestring(self) -> str:
        str_list = [
            '!'+'='*49,
            'BODY',
            '%s' % self.name,
            '#Nbody\tBspace',
            '%d %f' % (self.Nbody, self.Bspace)
            ]

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

        if self.X1 is not None and self.X2 is not None:
            str_list.extend([
                '#',
                'BFILE\t %f %f' % (self.X1, self.X2),
                self.Bfile
            ])
        else:
            str_list.extend([
                '#',
                'BFILE',
                self.Bfile
            ])
        return '\n'.join(str_list)


    @classmethod
    def _fromstring(cls,string):
        lines = string.split('\n')
        name = lines.pop(0)
        temp = lines.pop(0).split()
        Nbody,Bspace = int(float(temp[0])),float(temp[1])

        Ydupl = None
        scale = None
        translate = None
        Bfile,X1,X2 = None,None,None

        while lines:
            line = lines.pop(0)
            keyw = line[0:4]
            if keyw == 'YDUP':
                line = lines.pop(0)
                temp = line.split()
                Ydupl = float(temp[0])
            elif keyw == 'SCAL':
                line = lines.pop(0)
                temp = line.split()
                scale = (float(temp[0]),float(temp[1]),float(temp[2]))
            elif keyw == 'TRAN':
                line = lines.pop(0)
                temp = line.split()
                translate = (float(temp[0]),float(temp[1]),float(temp[2]))
            elif keyw == 'BFIL':
                line = lines.pop(0)
                temp = line.split()
                if len(temp) >= 3:
                    Bfile,X1,X2 = temp[0],float(temp[1]),float(temp[2])
                else:
                    Bfile = temp[0]

        return cls(name,Nbody,Bspace,translate,Bfile,Ydupl,scale,X1,X2)
        
@dataclass
class configuration:
    #attr   dtype   default
    header: header
    surfaces: list[surface]
    bodies: list[body]

    @property
    def configfilestring(self) -> str:
        str_list = [self.header.configfilestring]
        for body in self.bodies:
            str_list.append(body.configfilestring)
        for surface in self.surfaces:
            str_list.append(surface.configfilestring)
        str_list.append('')
        return '\n'.join(str_list)

    def tofile(self,fname):
        with open(fname,'w') as file:
            file.write(self.configfilestring)
        return

    @classmethod
    def fromfile(cls,fname):
        comments = compile('[#|!].*')
        spaces = compile(' +')
        first_leading_space = compile('^\s*\n')
        leading_space = compile('\n\s+')
        empty = compile('\n\s*\n')
        SURF = compile('SURF.*?\n|BODY.*?\n')

        with open(fname,'r') as file:
            configfilestring = file.read()

        configfilestring = sub(spaces,' ',configfilestring)
        configfilestring = sub(comments,'\n',configfilestring)
        configfilestring = sub(empty,'\n',configfilestring)
        configfilestring = sub(first_leading_space,'',configfilestring)
        configfilestring = sub(leading_space,'\n',configfilestring)

        configfilesegments = split(SURF,configfilestring)
        segmenttypes = findall(SURF,configfilestring)
        avlheader = header._fromstring(configfilesegments.pop(0))
        surfaces = list()
        bodies = list()
        for type,segment in zip(segmenttypes,configfilesegments):
            if type[0:4] == 'SURF':
                surfaces.append(surface._fromstring(segment))
            if type[0:4] == 'BODY':
                bodies.append(body._fromstring(segment))
        return cls(avlheader,surfaces,bodies)



if __name__ == "__main__":
    from os import listdir
    from os.path import isfile, join

    avlconf = configuration.fromfile('runs/navion_fancy.avl')

    mainwing = avlconf.surfaces[0]

    mainwing.splitspan(6.0,0.0)

    avlconf.bodies = []
    avlconf.tofile('test.avl')


    #avlconf = configuration.fromfile('runs/navion_fancy.avl')
    #avlconf.tofile('test.avl')
    #files = [f for f in listdir('runs') if isfile(join('runs', f))]
    #for file in files:
    #    if file[-4:] == '.avl':
    #        print(file)
    #        avlconf = configuration.fromfile('runs/'+file)



    # TODO Claf is a section parameter and not an airfoil parameter
    #add run function

    

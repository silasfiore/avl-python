import avl
import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass,field,make_dataclass

DEG2RAD = math.pi/180.0


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

def vsp2avl(fname):
    xmlroot = ET.parse(fname).getroot()
    Vehicle = recparse(xmlroot.find("./Vehicle"))
    surfaces = list()

    # ============================================================
    #   GENERATE HEADER
    # ============================================================

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
        elif Geom.GeomBase.TypeName == 'Blank':
            if Geom.ParmContainer.Name == 'CG':
                XYZref = (  Geom.ParmContainer.XForm.X_Location,\
                            Geom.ParmContainer.XForm.Y_Location,\
                            Geom.ParmContainer.XForm.Z_Location)
    
    header = avl.header(name,Sref,Cref,Bref,XYZref)

    
    for Geom in Vehicle.Geom:
        # ============================================================
        #   GENERATE SURFACES
        # ============================================================
        if Geom.GeomBase.TypeName == 'Wing':
            #-------------------------------------------------
            # find Ydupl for mirrored wings
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
            #-------------------------------------------------
            # read other surface parameters
            surfname = Geom.ParmContainer.Name

            Nchord = 9
            Cspace = 1.0
            Nspan = 22
            Sspace = 2.0

            sections = list()

            translate = [   Geom.ParmContainer.XForm.X_Location,\
                            Geom.ParmContainer.XForm.Y_Location,\
                            Geom.ParmContainer.XForm.Z_Location]

            Xrot = Geom.ParmContainer.XForm.X_Rotation
            cx = math.cos(Xrot*DEG2RAD)
            sx = math.sin(Xrot*DEG2RAD)
            Yrot = Geom.ParmContainer.XForm.Y_Rotation
            cy = math.cos(Yrot*DEG2RAD)
            sy = math.sin(Yrot*DEG2RAD)

            reldihedralflag = bool(int(Geom.ParmContainer.WingGeom.RelativeDihedralFlag))
            reltwistflag = bool(int(Geom.ParmContainer.WingGeom.RelativeTwistFlag))

            #-------------------------------------------------
            # read Xsections
            for idx,XSec in enumerate(Geom.WingGeom.XSecSurf.XSec):
                XSecParams = XSec.ParmContainer.XSec
                if idx == 0:
                    Xle = 0.0
                    Yle = 0.0
                    Zle = 0.0
                    ainc = XSecParams.Twist

                    #correct shift of root le due to inclination of wing
                    dx = (1-math.cos(ainc*DEG2RAD))*XSecParams.Twist_Location*XSecParams.Tip_Chord
                    dy = 0
                    dz = (math.sin(ainc*DEG2RAD))*XSecParams.Twist_Location*XSecParams.Tip_Chord

                    translate[0] +=   cy*dx       + 0*dy      + sy*dz
                    translate[1] +=   sx*sy*dx    + cx*dy     - sx*cy*dz
                    translate[2] +=  -cx*sy*dx    + sx*dy     + cx*cy*dz

                else:
                    #leading edge sweep
                    if bool(int(XSecParams.Sweep_Location)):
                        sweeple = XSecParams.Sweep_Location
                    else:
                        sweeple = math.atan(math.tan(XSecParams.Sweep*DEG2RAD) + ( 2.0*XSecParams.Sweep_Location)*(1-XSecParams.Taper)/(XSecParams.Aspect * ( 1 + XSecParams.Taper)))

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

                    dx = math.tan(sweeple)*XSecParams.Span
                    dy = XSecParams.Span*math.cos(absoluteDihedral*DEG2RAD)
                    dz = XSecParams.Span*math.sin(absoluteDihedral*DEG2RAD)

                    Xle +=  cy*dx       + 0*dy      + sy*dz
                    Yle +=  sx*sy*dx    + cx*dy     - sx*cy*dz
                    Zle += -cx*sy*dx    + sx*dy     + cx*cy*dz

                    ainc += relativeTwist

                chord = XSec.ParmContainer.XSec.Tip_Chord

                if XSec.XSec.XSecCurve.XSecCurve.Type == '7':
                    #NACA PROFILE
                    tovc = XSec.XSec.XSecCurve.ParmContainer.XSecCurve.ThickChord
                    camber = XSec.XSec.XSecCurve.ParmContainer.XSecCurve.Camber
                    camberloc = XSec.XSec.XSecCurve.ParmContainer.XSecCurve.CamberLoc
                    type = avl.aftypes.NACA
                    X1 = 0.0
                    X2 = 1.0
                    NACA = '%d%d%02d' % (round(camber*100),round(camberloc*10),round(tovc*100))

                    af = avl.airfoil(type,X1,X2,NACA,None,None,None,None)
                else:
                    af = None

                sections.append(avl.section((Xle,Yle,Zle),chord,ainc,None,None,af,None,None))
            
            surfaces.append(avl.surface(surfname,Nchord,Cspace,Nspan,Sspace,sections,tuple(translate),None,Ydupl))

        # ============================================================
        #   GENERATE BODIES
        # ============================================================            
            




    config = avl.configuration(header,surfaces,[])

    return config


if __name__ == "__main__":
    config = vsp2avl('glider.vsp3')

    for section in config.surfaces[0].sections:
        print(section)
    config.tofile('misc.avl')
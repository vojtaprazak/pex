def get_surface_atoms_with_mod9v3(struct, cutoff):
    struct.writeToPDB('monkey.pdb')
    directory = os.getcwd()

    lines = """# Example for: model.write_data()

from modeller import *
from modeller.scripts import complete_pdb
log.verbose()

env = environ()
env.io.atom_files_directory = '"""+directory+"""'
env.io.hetatm = False
env.io.water = False

env.libs.topology.read(file='$(LIB)/top_heav.lib')
mdl = model(env)
mdl.read(file='monkey.pdb')

myedat = energy_data()
myedat.radii_factor = 1.0 # The default is 0.82 (for soft-sphere restraints)
mdl.write_data(file='monkey', edat=myedat, output='PSA')"""

    f = file('write_data2.py', 'w')
    f.write(lines)
    f.close()
    system('mod9v3 write_data2.py')
    a = Structure('monkey.sol', surface = True)
    b = []
    for x in a.atomList:
        if x.access > cutoff:
            b.append(x)
    return Structure(b)
        

def makeVector(atom):
    return Vector(atom.getX(), atom.getY(), atom.getZ())

def translateToAtom(vector, atom):
    newAtom = Atom(atom.allInfo[:])
    newAtom.setX(vector.x)
    newAtom.setY(vector.y)
    newAtom.setZ(vector.z)
    return newAtom

def getSimiliarVectors(vectorList, vector, cutoff):
    similiarVectors = []
    for v in vectorList:
        score = abs(v.x-vector.x) + abs(v.y-vector.y) + abs(v.z-vector.z)
        if(score < cutoff):
            similiarVectors.append(v)
    return similiarVectors
    
    
def getPhiAngles(angleList):
    phis = []
    for x in range(len(angleList)):
        if(x%3==2):
            phis.append(angleList[x])
    return phis
    
def getPsiAngles(angleList):
    psis = []
    for x in range(len(angleList)):
        if(x%3==0):
            psis.append(angleList[x])
    return psis

def getOmegaAngles(angleList):
    omegas = []
    for x in range(len(angleList)):
        if(x%3==1):
            omegas.append(angleList[x])
    return omegas
        
def rebuildAngleList(phis, psis, omegas):
    angleList = []
    for x in range(len(psis)):
        angleList.append(psis[x])
        angleList.append(omegas[x])
        angleList.append(phis[x])
    return angleList

def getBondAngles(aList):
    bondAngles = []
    for x in range(len(aList)-2):
        A = makeVector(aList[x])
        B = makeVector(aList[x+1])
        C = makeVector(aList[x+2])
        AB = B.minus(A)
        BC = C.minus(B)
        bondAngles.append(acos((AB.dot(BC))/(AB.mod()*BC.mod())))
    return bondAngles

def placeAtom(A, B, C, l, ang, torsion):
    A = makeVector(A)
    B = makeVector(B)
    C = makeVector(C)

    AB = B.minus(A)
    BC = C.minus(B)
    
    n = AB.cross(BC).unit()
    BC = BC.unit()
    nxBC = n.cross(BC)
    
    dummyD = Vector(l*cos(ang), l*cos(torsion)*sin(ang), l*sin(torsion)*sin(ang))
    
    newX = dummyD.x*BC.x + dummyD.y*nxBC.x + dummyD.z*n.x + C.x
    newY = dummyD.x*BC.y + dummyD.y*nxBC.y + dummyD.z*n.y + C.y
    newZ = dummyD.x*BC.z + dummyD.y*nxBC.z + dummyD.z*n.z + C.z
    
    return Vector(newX, newY, newZ)

def rebuildBackbone(atomList, torsionAngleList):
    newAtomList = []
    for x in range(3):
        newAtomList.append(atomList[x].copy())
    for x in range(len(torsionAngleList)):
        lA = getLenAng(atomList[x+1].copy(), atomList[x+2].copy(), atomList[x+3].copy())
        v = placeAtom(newAtomList[x],newAtomList[x+1],newAtomList[x+2], lA[0], lA[1], torsionAngleList[x])
        #print v
        newAtomList.append(translateToAtom(v, atomList[x+3].copy()))
    return newAtomList
    
def rebuildStructure(fullAtomList, torsionAngleList):
    newAtomList = []
    oldBackbone = getBackbone(fullAtomList)
    oldTorsionAngs = torsionAngles(vectorise(oldBackbone))
    newBackbone = rebuildBackbone(oldBackbone, torsionAngleList)
    oldSideChains = getSideChains(fullAtomList)
    phiAngs = getPhiAngs(torsionAngleList)
    for x in range(0, len(newBackbone), 3):
        newAtomList.append(newBackbone[x])
        newAtomList.append(newBackbone[x+1])
        newAtomList.append(newBackbone[x+2])
        resNo = newBackbone[x].getRes()
        for y in oldSideChains:
            if y.getRes() == resNo:
                o = rotationFootpoint(newBackbone[0], newBackbone[1], y)
                
        
def getSideChains(atomList):
    sideChains = []
    for x in atomList:
        if(x.allInfo[2] == 'CA' or x.allInfo[2] == 'N' or x.allInfo[2] == 'C' or x.allInfo[2] == 'O'):
            pass
        else:
            sideChains.append(x.copy())
    return sideChains
        

def topCCD(firstAtom, secondAtom, mobile, anchor):
    fst = makeVector(firstAtom)
    snd = makeVector(secondAtom)
    m = makeVector(mobile)
    a = makeVector(anchor)
    
    o = rotationFootpoint(fst, snd, m)
    f = a.minus(o) 
    axis = snd.minus(fst)
    r = m.minus(o)
    s = (r.unit()).cross(axis.unit())

    return -(f.dot(s))*r.mod()

def bottomCCD(firstAtom, secondAtom, mobile, anchor):
    fst = makeVector(firstAtom)
    snd = makeVector(secondAtom)
    m = makeVector(mobile)
    a = makeVector(anchor)
    o = rotationFootpoint(fst, snd, m)
    f = a.minus(o) 
    r = m.minus(o)
    return f.dot(r.unit())*r.mod()

def getCCDAngle(firstAtom, secondAtom, mobiles, anchors):
    top = 0
    bottom = 0
    sndDer1 = 0
    sndDer2 = 0
    fst = makeVector(firstAtom)
    snd = makeVector(secondAtom)
    axis = snd.minus(fst)
    for x in range(len(mobiles)):
        m = makeVector(mobiles[x])
        a = makeVector(anchors[x])
        o = rotationFootpoint(fst, snd, m)
        f = a.minus(o) 
        r = m.minus(o)
        s = (r.unit()).cross(axis.unit())
        top += -(f.dot(s))*r.mod()
        bottom += f.dot(r.unit())*r.mod()
        sndDer1 += 2*r.mod()*(f.unit().dot(r.unit()))
        sndDer2 += 2*r.mod()*(f.unit().dot(s))
    angle = atan(top/bottom)
    sndDer = sndDer1*cos(angle) + sndDer2*sin(angle)
    if sndDer > 0:
        return angle
    else:
        return angle-sign(angle)*pi
        
def sign(number):
    if(number < 0):
        return -1
    elif(number > 0):
        return 1
    else:
        return 0        
        
def runCCD(iterations, beforeBreak, breakSection, afterBreak):
    anchors = []
    mobiles = []
    anchors.append(afterBreak.getAtom(0))
    anchors.append(afterBreak.getAtom(1))
    anchors.append(afterBreak.getAtom(2))
    ccdList = breakSection.getAtomList()
    mobiles.append(breakSection.getAtom(-3))
    mobiles.append(breakSection.getAtom(-2))
    mobiles.append(breakSection.getAtom(-1))
    for x in range(iterations):
        for x in range(len(ccdList)-4):
            if(ccdList[x].getName() == 'CA'):
                pass
            else:
                mobiles[0] = ccdList[-3]
                mobiles[1] = ccdList[-2]
                mobiles[2] = ccdList[-1]
                tAngles = (Structure(ccdList)).getTorsionAngles()
                firstAtom = ccdList[x+1]
                secondAtom = ccdList[x+2]
                newAngle = getCCDAngle(firstAtom, secondAtom, mobiles, anchors)
                tAngles[x] += newAngle
                ccdList = rebuildBackbone(ccdList, tAngles)
                distance = (mobiles[0].minus(anchors[0])**2 + mobiles[1].minus(anchors[1])**2 + mobiles[2].minus(anchors[2])**2)**0.5
                print(distance)
                if(distance < 0.08):
                    return Structure(beforeBreak.getAtomList() + ccdList + afterBreak.getAtomList()[3:])
    return Structure(beforeBreak.getAtomList() + ccdList + afterBreak.getAtomList()[3:])
                
def getLenAng(fstAtom, sndAtom, thdAtom):
    B = makeVector(fstAtom)
    C = makeVector(sndAtom)
    D= makeVector(thdAtom)
    l = D.minus(C).mod()
    ang = (D.minus(C)).arg(C.minus(B))
    return (l, ang)

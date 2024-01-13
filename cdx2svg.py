import struct
import math
import re

pt={1:"H",2:"He",3:"Li",4:"Be",5:"B",6:"C",7:"N",8:"O",9:"F",10:"Ne",11:"Na"}
pal=[]
nre=re.compile("((<sub>)?(\d+)(</sub>)?)")
nrc=re.compile("((<sub>)?([\+-])(</sub>)?)")

class obj:
  def __repr__(self):
   return vars(self).__str__()+"\n"

def formula_format(s): # add <sub> to numbers, <sup> to +/-, if not present already   
 return re.sub(nre,"<sub>\g<3></sub>",s)
   
def getNumber(d): # get number, figuring the length
 if len(d)==1: return struct.unpack('=B',d)[0]
 if len(d)==2: return struct.unpack('=H',d)[0]
 if len(d)==4: return struct.unpack('=L',d)[0]
 assert 1, "Len is not 1,2, or 4"

def getColor(n):
 c=pal[n-2]
 return "#{0:02x}{1:02x}{2:02x}".format(c[0],c[1],c[2])
 
def export2svg(objs,fname="a.cdml"):
 with open(fname,"w") as f:
  f.write('<?xml version="1.0" ?><cdml version="0.15" xmlns="http://www.freesoftware.fsf.org/bkchem/cdml">\n')
  f.write(' <info>\n  <author_program version="0.13.0">BKchem</author_program>\n </info>\n')
  f.write(' <paper crop_margin="10" crop_svg="1" orientation="portrait" replace_minus="0" type="A4" use_real_minus="0"/>\n')
  f.write(' <viewport viewport="0.000000 0.000000 640.000000 480.000000"/>\n')
  f.write(' <standard area_color="" font_family="helvetica" font_size="12" line_color="#000" line_width="1px" paper_crop_margin="10" paper_crop_svg="0" paper_orientation="portrait" paper_type="A4">\n')
  f.write('  <bond double-ratio="0.75" length="0.7cm" wedge-width="5px" width="6px"/>\n')
  f.write('   <arrow length="1.6cm"/>\n')
  f.write('   <atom show_hydrogens="0"/>\n')
  f.write(' </standard>\n')


  # process data: make list of atoms and bond, calculate coordinates

  atoms={}
  bonds=[]
  graphx=[]
  txts=[]

  scale=1.2
  
  for ob in objs:
   if not hasattr(ob,"t"): 
    print "No t",ob
    continue
   
   if ob.t=="m": pass
   elif ob.t=="a":    
    catom=ob
    catom.label="C"
    catom.x=round(ob.XY[1]/1857710.0*scale,4)
    catom.y=round(ob.XY[0]/1857710.0*scale,4)
    atoms[ob.id]=catom
    if hasattr(ob,"tobj"): #text
     pass
   elif ob.t=="b":
    bonds.append(ob)
   elif ob.t=="t": # text
    if ob.id:
     txts.append(ob)
     ob.x=round(ob.XY[1]/1857710.0*scale,4)
     ob.y=round(ob.XY[0]/1857710.0*scale,4)
   elif ob.t=="g": # graphic, 
    graphx.append(ob)
    if ob.grtype==1 or ob.grtype==3: #line or rectangle
     if not hasattr(ob,"fcolor"): ob.fcolor=0
     ob.y1=round(ob.bbox[0]/1857710.0*scale,4)
     ob.x1=round(ob.bbox[1]/1857710.0*scale,4)
     ob.y2=round(ob.bbox[2]/1857710.0*scale,4)
     ob.x2=round(ob.bbox[3]/1857710.0*scale,4)
    elif ob.grtype==2: # arc
     if not hasattr(ob,"fcolor"): ob.fcolor=0
     ob.x1=round(ob.bbox[1]/1857710.0*scale,4)
     ob.y1=round(ob.bbox[0]/1857710.0*scale,4)
     cx=round(ob.bbox[3]/1857710.0*scale,4)
     cy=round(ob.bbox[2]/1857710.0*scale,4)
     dx=ob.x1-cx
     dy=ob.y1-cy
     ang=ob.arcsize/1800*math.pi
     print "ang", ang
#     ob.x3=cx-dx
#     ob.y3=cy-dy
#     ob.x2=cx
#     ob.y2=cy

     ob.x3=round(cx+math.cos(ang)*dx-math.sin(ang)*dy,4)
     ob.y3=round(cy+math.sin(ang)*dx+math.cos(ang)*dy,4)
     ob.x2=round(cx+math.cos(ang/2)*dx-math.sin(ang/2)*dy,4)
     ob.y2=round(cy+math.sin(ang/2)*dx+math.cos(ang/2)*dy,4)
    elif ob.grtype==4: pass  # oval
    elif ob.grtype==5: pass  # orbital
    elif ob.grtype==6: pass  # bracket
    elif ob.grtype==7: pass # symbol

  mols=[]

  # go though all bonds to find connected fragments
  for b in bonds:   
   ms=[]
   for m in mols:
    if b.fromatom in m.atomids or b.fromatom in m.atomids:
     m.bonds.append(b)
     m.atomids.add(b.fromatom)
     m.atomids.add(b.toatom)
     ms.append(m)
   
   if len(ms)==0: # no connections found, start new molecule
    m=obj()
    m.bonds=[b]
    m.atomids=set((b.fromatom,b.toatom))
    mols.append(m)
   elif len(ms)>1: # connections to multiple molecules found, merge them
    m1=ms.pop()
    m1.append(b)
    m.atomids=set((b.fromatom,b.toatom))
    for m in ms:
     m1.bonds.extend(m.bonds)
     m1.atomids.update(m.atomids)
     mols.remove(m)
  
  print "mols",len(mols)  
  # check for orphan atoms, add to extra molecule
  aa=set()
  for m in mols:
   aa.update(m.atomids)
  oa=aa.difference(set(atoms.keys())) 
  
  if oa:
   m=obj()
   m.bonds=[]
   m.atomids=oa
   mols.append(m)

  # output molecules
  n=0
  for m in mols:
   n+=1
   f.write(' <molecule id="molecule%s" name="">\n' % n)
   for a_id in m.atomids:
    ob=atoms[a_id]
    e=""
    e1=""
    for a,v in ob.__dict__.items():
     if a=="anum": ob.label=pt[v]
     elif a=="Hs" and v: e+=' hydrogens="on"'
     elif a=="charge": e+=' charge="%s"' % v
     elif a=="fcolor": e1+='<font color="%s" family="helvetica" size="12"/>' % getColor(v)

   #<mark auto="0" draw_circle="yes" size="10.0" type="plus" x="1.551cm" y="3.072cm"/>
    if hasattr(ob,"tobj"):
     f.write('  <text id="atom%s" pos="center-first">%s<ftext>%s</ftext><point x="%scm" y="%scm"/></text>\n' % (ob.id,e1,formula_format(ob.tobj.text),ob.x,ob.y))
    else:
     f.write('  <atom id="atom%s" name="%s"%s>%s<point x="%scm" y="%scm"/></atom>\n' % (ob.id,ob.label,e,e1,ob.x,ob.y)) # valency="4"
    
   for ob in m.bonds:
    l="n1"
    e=""
    if hasattr(ob,"border"):
      if ob.border==2: l="n2"
    if hasattr(ob,"display"):
      if ob.display==3: l="h1"
      if ob.display==6: l="w1"       
      elif ob.display==4 or ob.display==7: 
       l="h1" if obj.display==4 else "w1"
       ob.fromatom, ob.toatom=ob.toatom,ob.fromatom

    if hasattr(ob,"fcolor"): e+=' color="%s"' % getColor(ob.fcolor)
       
    f.write('  <bond double_ratio="0.75" id="bond%s" line_width="1.0" start="atom%s" end="atom%s" type="%s" %s/>\n' % (ob.id,ob.fromatom,ob.toatom,l,e))
   f.write(' </molecule>\n')  

  for ob in txts:
    f.write('<text id="text%s"><point x="%scm" y="%scm"/><ftext>%s</ftext></text>' % (ob.id,ob.x,ob.y,ob.text))
  
  
  for ob in graphx:
   if ob.grtype==1: #line
    if hasattr(ob,"arrowtype"):
     start="no"
     end="no"
#       0 - no head, 1 - half, 2 - reg, 32 - retro?, 16 - hollow, 8 - eq?, 4 - res, 64 - nogo, 128 - dipole       
    if ob.arrowtype==1: pass # half
    elif ob.arrowtype==2 or ob.arrowtype==4 or ob.arrowtype==8: # regular
     start="no" if ob.arrowtype==4 else "yes" 
     end="no" if ob.arrowtype==2 or ob.arrowtype==4 else "yes"
    f.write('  <arrow color="%s" id="arrow%s" shape="(8, 10, 3)" spline="no" start="%s" end="%s" type="normal" width="1.0">\n  <point x="%scm" y="%scm"/>\n  <point x="%scm" y="%scm"/>\n  </arrow>\n' % (getColor(ob.fcolor),ob.id,start,end,ob.x1,ob.y1,ob.x2,ob.y2))
   elif ob.grtype==2: # arc
    f.write('  <arrow color="%s" id="arrow%s" shape="(8, 10, 3)" spline="yes" start="%s" end="%s" type="normal" width="1.0">\n  <point x="%scm" y="%scm"/>\n  <point x="%scm" y="%scm"/>\n  <point x="%scm" y="%scm"/>\n  </arrow>\n' % (getColor(ob.fcolor),ob.id,"yes","no",ob.x1,ob.y1,ob.x2,ob.y2,ob.x3,ob.y3))
   elif ob.grtype==3: # rectangle
     f.write('<rect area_color="" line_color="#000" width="1.0" x1="%scm" x2="%scm" y1="%scm" y2="%scm"/>' % (ob.bbxo[0],ob.bbxo[1],ob.bbxo[2],ob.bbxo[3]))
   elif ob.grtype==4: pass  # oval
   elif ob.grtype==5: pass  # orbital
   elif ob.grtype==6: pass  # bracket
   elif ob.grtype==7: pass # symbol
     
  f.write('</cdml>')

def readcdx(f):
 if f.read(8)!="VjCD0100": return
 f.read(4)#4,3,2,1
 f.read(16)# 0s
 objs=[]
 objstack=[]
 curobj=None
 prevobj=None
 
 while True:
  tag=f.read(2)
  if tag is None: break
  tag=struct.unpack('=H',tag)[0] # tag
  print "tag","0x%x" % tag
  if tag & 0x8000: # object
   id=struct.unpack('=l',f.read(4))[0] #id
   print "id",id
   if curobj:
    objstack.append(curobj) 
   prevobj=curobj
   curobj=obj()
   curobj.id=id
   objs.append(curobj)
   
   if tag==0x8003: # molecule
    curobj.t="m"
   elif tag==0x8004: # atom
    curobj.t="a"
   elif tag==0x8005: # bond
    curobj.t="b"
   elif tag==0x8006: # text
    curobj.t="t"
    if prevobj: prevobj.tobj=curobj
   elif tag==0x8007: # graphic
    curobj.t="g"
   elif tag==0x8008: # curve
    curobj.t="c"
   elif tag==0x8021: # graphic
    curobj.t="G"
   elif tag==0x8027: # arrow
    curobj.t="A"
  else: #property
   if tag==0:
    try:
     print "End object",curobj.t,curobj.id
     curobj=objstack.pop()
    except:
     print objs
     export2svg(objs)
     exit(0)
    continue
   
   size=struct.unpack('=H',f.read(2))[0]
   if size==0xFFFF: size=struct.unpack('=L',f.read(4)[0])
   print "size",size
   data=f.read(size)

   if tag==0x0200: # X,Y
    curobj.XY=struct.unpack('=LL',data)
   elif tag==0x0204: # 
    n=struct.unpack('=LLLL',data)
    if curobj: curobj.bbox=n
   elif tag==0x0300: # palette
    n=struct.unpack('=H',data[:2])[0]
#    pal=[]
    for i in range(n):
     c=struct.unpack('=HHH',data[2+i*6:8+i*6])
     pal.append((c[0]/256,c[1]/256,c[2]/256))
    print pal
   elif tag==0x0301: # fcolor
    if curobj: curobj.fcolor=struct.unpack('=H',data[:2])[0]    
   elif tag==0x0302: # bcolor
    if curobj: curobj.bcolor=struct.unpack('=H',data[:2])[0]    
   elif tag==0x0402: # text
    curobj.anum=struct.unpack('=H',data)[0]
   elif tag==0x0421: # charge
    curobj.charge=getNumber(data) # struct.unpack('=b',data[0])[0]
   elif tag==0x042B: # Hs
    curobj.Hs=struct.unpack('=H',data)[0]
   elif tag==0x042E: # stereo
    curobj.stereo=struct.unpack('=H',data)[0]
   elif tag==0x0604:  # from atom
    curobj.fromatom=struct.unpack('=L',data)[0]   
   elif tag==0x0605:  # to atom
    curobj.toatom=struct.unpack('=L',data)[0]
   elif tag==0x0600: # bond order
    curobj.border=struct.unpack('=H',data)[0]
   elif tag==0x0601:  # bond display 
    curobj.display=struct.unpack('=H',data)[0]   # 1 dash, 2 hash, 3,4-dash b/e, 6,7 - wedge b/e
   elif tag==0x0602: pass  # 2nd bond display
   elif tag==0x0603: # dbposition
    curobj.dbpos=struct.unpack('=H',data)[0]   
   elif tag==0x060A: #stereo
    curobj.stereo=struct.unpack('=B',data)[0]
   elif tag==0x060B: #bond list
    curobj.blist=struct.unpack('=LLLL',data)
   elif tag==0x0700: # text
    ns=struct.unpack('=H',data[:2])[0]
    curobj.text=data[2+ns*10:]
#    for i in range(ns):
#     st=struct.unpack('=HHHHH',data[2+i*10:2+i*10+10])
   elif tag==0x0701: # h justification
    curobj.hjust=struct.unpack('=B',data)[0]
   elif tag==0x0702: # line height
    curobj.lheight=struct.unpack('=H',data)[0]
   elif tag==0x0705: # node alignment
    curobj.nalign=struct.unpack('=B',data)[0]
   elif tag==0x0706: # text line height
    curobj.tlheight=struct.unpack('=H',data)[0]
   elif tag==0x0A00: # gr type 
    curobj.grtype=struct.unpack('=B',data[0])[0]    #  1-line,2-arc,3 -rectangle,  4-oval,5-orbital, 6 -bracket 7 - symbol
   elif tag==0x0A02: # arrow type
    curobj.arrowtype=struct.unpack('=B',data[0])[0] #  0 - no head, 1 - half, 2 - reg, 32 - retro?, 16 - hollow, 8 - eq?, 4 - res, 64 - nogo, 128 - dipole
   elif tag==0x0A20: # arrow head size
    curobj.arrowhsize=struct.unpack('=H',data)[0]
   elif tag==0x0A21: # arc size
    curobj.arcsize=struct.unpack('=H',data)[0]
   elif tag==0x081A:  pass # def font fam
   elif tag==0x081C:  pass # def fomt size
   elif tag==0x081E:  pass # def fomt face
   

f=open("a.cdx","r")
readcdx(f)

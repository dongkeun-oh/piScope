import wx, os, re, string, logging, shutil
from ifigure.mto.py_contents import Efitgfile
from numpy import *
from collections import OrderedDict

efitnodes =['\\bcentr','\\bdry','\\cpasma','\\epoten','\\ffprim',
            '\\fpol','\\gtime','\\lim','\\limitr','\\mh','\\mw','\\nbdry',
            '\\pprime','\\pres','\\psin','\\psirz','\\qpsi','\\r','\\rgrid1',
            '\\rhovn','\\rmaxis','\\rzero','\\ssibry','\\ssimag','\\xdim','\\z',
            '\\zdim','\\zmaxis','\\zmid']

def string_split_ig(s, sep):
   t=re.split(sep,s)
   id=range(len(t))
   id.reverse()
   for i in id:
      if t[i]=='': del t[i]
   return t

def eval_form2020(s):
   #form2020= '(5e16.9)'
   d=16
   ret=[]
   for i in range(5) :
#     print s[d*i:d*(i+1)],s[d*i:d*(i+1)].strip().__repr__()
     if s[d*i:d*(i+1)].strip() != '':
        ret.append(float(s[d*i:d*(i+1)]))     
   return ret

def readandstrip(f):
   line = ''
   while line == '':
     line=f.readline()
     line = line.rstrip("\r\n")
   return line

def make_form2020(array):
   line0=[]
   line = ''
   i=0
   while i < len(array):
      line=line+'{0: 13.9e}'.format(array[i])
      i = i+1
      if (i % 5) == 0:
         line0.append(line)
         line = ''

   if len(line) != 0:
      line0.append(line)

   return line0

def load_array_form2020(f, size):
   k=0
   arr=[]
   while (k < size):
     arr=arr+eval_form2020(readandstrip(f))
     k = k+5
   return array(arr)

def load_matrix_form2020(f, size1, size2):
   k=0
   arr=[]
   while k < size1*size2:
     arr=arr+eval_form2020(readandstrip(f))
     k = k+5   
   return array(arr).reshape(size1, size2)

def load_file(file=None):
   def is_ascii(s):
       try:
           s.encode().decode('ascii')
       except:
           return False
       return True
#       return all([ord(c) < 128 for c in s])
#   file='/home/shiraiwa/g1101019014.01340'
   f = open(file, 'r')
   
   line=f.readline()
   line = line.rstrip("\r\n")
   tmp=string_split_ig(line[49:], ' |,')
   header=line[:48]
   idum=int(tmp[0])
   mw=int(tmp[1])
   mh=int(tmp[2])

   xdim,zdim,rzero,rgrid1,zmid=eval_form2020(readandstrip(f))
   rmaxis,zmaxis,ssimag,ssibdry,bcentr=eval_form2020(readandstrip(f))
   cpasma,ssimag,xdum,rmaxis,xdum=eval_form2020(readandstrip(f))
   zmaxis,xdum,ssibry,xdum,xdum=eval_form2020(readandstrip(f))

   fpol=load_array_form2020(f, mw)
   pres=load_array_form2020(f, mw)
   ffprim=load_array_form2020(f, mw)
   pprime=load_array_form2020(f, mw)
   
   psirz=load_matrix_form2020(f, mh, mw)
   qpsi=load_array_form2020(f, mw)

   nbbbs,limitr=string_split_ig(readandstrip(f), ' |,')
   #print nbbbs, limitr
   nbbbs=int(nbbbs) 
   limitr=int(limitr)
   

   if nbbbs is not 0:
      rzbbbs = load_matrix_form2020(f, nbbbs,2)

   xylim=load_matrix_form2020(f, limitr, 2)

   rgrid = rgrid1 + xdim*arange(mw)/(mw-1.)
   zgrid = zmid-0.5*zdim + zdim*arange(mh)/(mh-1.)

   nm = Efitgfile()

   if cpasma > 0: 
     sss = -1
     ssimag = ssimag*sss
     ssibdry = ssibdry*sss
     ssibry = ssibry*sss
     psirz=psirz*sss

   nm["table"]={}
   val = nm["table"]
   val["header"] = header
   val["idum"]=idum
   val["mw"]=mw
   val["mh"]=mh
   val["xdim"]=xdim
   val["zdim"]=zdim
   val["rzero"]=rzero
   val["rgrid1"]=rgrid1
   val["zmid"]=zmid
   val["rmaxis"] = rmaxis
   val["zmaxis"] = zmaxis
   val["ssimag"] = ssimag
   val["ssibdry"] = ssibdry
   val["bcentr"] = bcentr
   val["cpasma"] = cpasma
   val["ssibry"] = ssibry

   val["rgrid"]=rgrid
   val["zgrid"]=zgrid
   val["psirz"]=psirz
   val["fpol"] =fpol
   val["pres"] =pres
   val["ffprim"]=ffprim
   val["pprime"]=pprime
   val["qpsi"]=qpsi
   val["nbbbs"]=nbbbs
   val["rbbbs"]=rzbbbs[:,0]
   val["zbbbs"]=rzbbbs[:,1]
   val["nlim"]=limitr
   val["xlim"]=xylim[:,0]
   val["ylim"]= xylim[:,1]

   # namelist section
   sec = None
   end_flag=False
   while 1:
      line0=f.readline()
      if not is_ascii(line0): continue
      if not line0: break
      line0 = line0.rstrip("\r\n")
      line0 = ' '.join(line0.split())
      if line0.endswith('/'):
         line = line0[:-1]
      else:
         line = line0

      if  line == '' : continue
      if line.startswith('&'):
         s=string_split_ig(line, ' |,')
         sec=s[0][1:]
         print( 'making new sec ', sec, line.__repr__(), s)
         nm[sec] = OrderedDict()
         if len(s) > 1:
            line = ' '.join(s[1:])
         else:
            continue
      if sec is None:
         continue #skip unitl first &*** starts
      sall=string_split_ig(line, ' |,')
 
      i=0
      while i < len(sall):
         if sall[i].find('=') != -1:
            if sall[i] == '=':    # '='
               sall[i-1]=sall[i-1]+'='
               del sall[i]
               continue
            if sall[i].startswith('='):  #'=value'
               sall[i-1]=sall[i-1]+'='
               sall[i]=sall[i][1:] 

               continue
            if sall[i].endswith('='):continue #'name='
            k=sall[i].split('=')
            sall[i]=k[0]+'='
            sall.insert(i+1, k[1])
         i=i+1

      for s in sall:
        if s.find('=') != -1:
          k=s.split('=')
          varname=k[0]
#          s[0]=k[1]
#          print s
#          if s[0] is '': del s[0]
          if debug != 0: print('create dict key', sec, varname)
          #print 'create dict key', sec, varname
          nm[sec][varname]=[]
          if s.endswith('/'):
             sec=None
          continue
      # for lines without 'xxxx = '
        #print s, sec, varname
        nm[sec][varname]=nm[sec][varname]+parseStr(s)
      if line0.endswith('/'):
           sec=None

   xx,yy = meshgrid(rgrid, zgrid) 
#   isPlasma = xx.copy()
   sss = len(rgrid)*len(zgrid)
   isPlasma = array([False]*sss).reshape(len(zgrid), len(rgrid))
   for m, xv in enumerate(rgrid):
       for n, yv in enumerate(zgrid):
           dx = rzbbbs[:,0] - xv
           dy = rzbbbs[:,1] - yv
           d1 = sqrt(dx[:-1]**2 + dy[:-1]**2)
           d2 = sqrt(dx[1:]**2 +  dy[1:]**2)

           d = (dx[:-1]*dy[1:] - dx[1:]*dy[:-1])/d1/d2
           d = d[abs(d) < 1.0]
           xxxx = sum(arcsin(d))
           isPlasma[n, m] = (abs(xxxx) > 3)
#              print isPlasma[n, m] > 3
#   val['isPlasma0'] = isPlasma
#   print 'here'
#   isPlasma =  abs(isPlasma) > 3

   xpsi = (psirz-ssimag)/(ssibdry-ssimag)
   for m, xv in enumerate(rgrid):
       for n, yv in enumerate(zgrid):
           if not isPlasma[n,m]:
               xpsi[n, m] = 1.0
   xxx = linspace(0, 1, fpol.shape[0])
   fp = interp(xpsi.flatten(), xxx, fpol)   
   fp = fp.reshape(xx.shape)
   fc = interp(xpsi.flatten(), xxx, ffprim)   
   fc = fc.reshape(xx.shape)
   pc = interp(xpsi.flatten(), xxx, pprime)
   pc = pc.reshape(xx.shape)   
   pr = interp(xpsi.flatten(), xxx, pres)
   pr = pr.reshape(xx.shape)
   qr = interp(xpsi.flatten(), xxx, qpsi)
   qr = qr.reshape(xx.shape)
   val["isPlasma"] = isPlasma

   pr[isPlasma != True] = 0.0
   pc[isPlasma != True] = 0.0
   fc[isPlasma != True] = 0.0
   val["pressrz"] = pr
   val["qrz"] = qr
   val["btrz"] = fp/xx

   dpsidz, dpsidr = gradient(psirz)
   brrz = -dpsidz/(zgrid[1]-zgrid[0])/xx
   bzrz = dpsidr/(rgrid[1]-rgrid[0])/xx
   val["brrz"] = brrz
   val["bzrz"] = bzrz

   mu0 = 4e-7*3.1415926535
   val["jtrz"] = -(xx*pc+fc/xx/mu0)/1e6  #1e6=(MA/m2)   
   
   
   k = (val["zmaxis"] - rgrid[0])/(rgrid[1] - rgrid[0])
   from scipy.interpolate import interp2d

   f = interp2d(rgrid, zgrid, psirz, kind = 'cubic')
   val['npsimid'] = array([(f(r, val["zmaxis"]) - ssimag)/(ssibry- ssimag) 
                              for r in rgrid]).flatten() 
   f1 = interp2d(rgrid, zgrid, val["btrz"], kind = 'cubic')
   f2 = interp2d(rgrid, zgrid, val["bzrz"], kind = 'cubic')
   val['gammamid'] =  array([arctan(f2(r, val["zmaxis"])
                                /f1(r, val["zmaxis"]))*180/3.1415926
                              for r in rgrid]).flatten()    
   val['bzmid'] = array([f2(r, val["zmaxis"])
                              for r in rgrid]).flatten() 
   val['btmid'] = array([f1(r, val["zmaxis"])
                              for r in rgrid]).flatten() 
   f = interp2d(rgrid, zgrid, val["jtrz"], kind = 'cubic')
   val['jtmid'] = array([f(r, val["zmaxis"])
                              for r in rgrid]).flatten() 
   f = interp2d(rgrid, zgrid, val["pressrz"], kind = 'cubic')
   val['pressmid'] = array([f(r, val["zmaxis"])
                              for r in rgrid]).flatten() 
   f = interp2d(rgrid, zgrid, val["qrz"], kind = 'cubic')
   val['qmid'] = array([f(r, val["zmaxis"])
                              for r in rgrid]).flatten() 
   return nm


def mdsdata2tree( mdsdata, tree, shot, idx):
   nm = Efitgfile()
   nm["table"]={}
   val = nm["table"]

   for node_str in efitnodes:
     name = node_str[1:]
     try:
        val[name] = mdsdata[name][idx]
     except:
        pass
        #print "mdsdata2tree : error in %s[%d]" %(name, idx)

   if mdsdata["cpasma"][idx] > 0. :  sgn = -1.
   else : sgn = 1.

   val["ssimag"] =sgn*val["ssimag"]
   val["ssibdry"] =sgn*val["ssibry"]  ##!!
   val["ssibry"] =sgn*val["ssibry"]
   val["psirz"] =sgn*val["psirz"]
   
   val["mw"] = int(val["mw"])
   val["mh"] = int(val["mh"])

   rgrid1 = val["rgrid1"]
   xdim = val["xdim"]
   zdim = val["zdim"]
   mw = val["mw"]
   mh = val["mh"]
   rgrid = rgrid1 + xdim*arange(mw)/(mw-1.)
   zgrid = val["zmid"]-0.5*zdim + zdim*arange(mh)/(mh-1.)
   val["rgrid"] = rgrid
   val["zgrid"] = zgrid
   nbbbs = val["nbdry"] ##!

   if nbbbs is not 0:
     rzbbbs = val["bdry"]
   xylim = mdsdata["lim"]
   val["nbbbs"] = int(nbbbs)
   val["rbbbs"] = rzbbbs[:int(nbbbs),0]
   val["zbbbs"] = rzbbbs[:int(nbbbs),1]
   val["nlim"] = int(val["limitr"])
   val["xlim"] = xylim[:,0]
   val["ylim"] = xylim[:,1]
 
# header
   import time
   gt=time.localtime(time.time())
   date='%02d/%02d/%4d'%(gt.tm_mon,gt.tm_mday,gt.tm_year)

   time_stamp = mdsdata["gtime"][idx]
#   header = " %-12s %10s #%6d %5dms         0"%(tree,date,shot,time_stamp)+str("%4i"%mw)+str("%4i"%mh)+"\n"
   header = " %-12s %10s #%6d %5dms         "%(tree,date,shot,time_stamp) 
   val["header"] = header
   val["idum"] = 0

   psirz = val["psirz"]
   ssimag = val["ssimag"]
   ssibdry = val["ssibdry"] 
   zmaxis = val["zmaxis"]

   xx, yy = meshgrid(rgrid, zgrid)
   sss = len(rgrid)*len(zgrid)
   isPlasma = array([False]*sss).reshape(len(zgrid), len(rgrid))
   for m, xv in enumerate(rgrid):
       for n, yv in enumerate(zgrid):
           dx = rzbbbs[:,0] - xv
           dy = rzbbbs[:,1] - yv
           d1 = sqrt(dx[:-1]**2 + dy[:-1]**2)
           d2 = sqrt(dx[1:]**2 +  dy[1:]**2)

           d = (dx[:-1]*dy[1:] - dx[1:]*dy[:-1])/d1/d2
           d = d[abs(d) < 1.0]
           xxxx = sum(arcsin(d))
           isPlasma[n, m] = (abs(xxxx) > 3)
   xpsi = (psirz-ssimag)/(ssibdry-ssimag)
   xxx = linspace(0, 1, val["fpol"].shape[0])
   fp = interp(xpsi.flatten(), xxx, val["fpol"])   
   fp = fp.reshape(xx.shape)
   fc = interp(xpsi.flatten(), xxx, val["ffprim"])   
   fc = fc.reshape(xx.shape)
   pc = interp(xpsi.flatten(), xxx, val["pprime"])
   pc = pc.reshape(xx.shape)   
   pr = interp(xpsi.flatten(), xxx, val["pres"])
   pr = pr.reshape(xx.shape)
   qr = interp(xpsi.flatten(), xxx, val["qpsi"])
   qr = qr.reshape(xx.shape)

   val["isPlasma"] = isPlasma
   pr[isPlasma != True] = 0.0
   pc[isPlasma != True] = 0.0
   fc[isPlasma != True] = 0.0
   val["pressrz"] = pr
   val["qrz"] = qr
   btrz = fp/xx
   val["btrz"] = btrz
    
   dpsidz, dpsidr = gradient(psirz)
   brrz = -dpsidz/(zgrid[1]-zgrid[0])/xx
   bzrz = dpsidr/(rgrid[1]-rgrid[0])/xx
   val["brrz"] = brrz
   val["bzrz"] = bzrz

   mu0 = 4e-7*3.1415926535
   val["jtrz"] = -(xx*pc+fc/xx/mu0)/1e6  #1e6=(MA/m2)   
   
   k = (zmaxis - rgrid[0])/(rgrid[1] - rgrid[0])

   from scipy.interpolate import interp2d

   f = interp2d(rgrid, zgrid, psirz, kind = 'cubic')
   val['npsimid'] = array([(f(r, zmaxis) - ssimag)/(ssibdry- ssimag) 
                              for r in rgrid]).flatten() 
   f1 = interp2d(rgrid, zgrid, btrz, kind = 'cubic')
   f2 = interp2d(rgrid, zgrid, bzrz, kind = 'cubic')
   val['gammamid'] =  array([arctan(f2(r, zmaxis)
                                /f1(r, zmaxis))*180/3.1415926
                              for r in rgrid]).flatten()    
   val['bzmid'] = array([f2(r, zmaxis)
                              for r in rgrid]).flatten() 
   val['btmid'] = array([f1(r, zmaxis)
                              for r in rgrid]).flatten() 
   f = interp2d(rgrid, zgrid, val["jtrz"], kind = 'cubic')
   val['jtmid'] = array([f(r, zmaxis)
                              for r in rgrid]).flatten() 
   f = interp2d(rgrid, zgrid, val["pressrz"], kind = 'cubic')
   val['pressmid'] = array([f(r, zmaxis)
                              for r in rgrid]).flatten() 
   f = interp2d(rgrid, zgrid, val["qrz"], kind = 'cubic')
   val['qmid'] = array([f(r, zmaxis)
                              for r in rgrid]).flatten() 
   return nm

def add_plotEQ( nmtable, viewer=None):
    x = nmtable["rgrid"]
    y = nmtable["zgrid"]
    z = nmtable["psirz"]
    xlim = nmtable["xlim"]
    ylim = nmtable["ylim"]

    rb = nmtable["rbbbs"]
    zb = nmtable["zbbbs"]
    fpol = nmtable["fpol"]
    ssimag = nmtable["ssimag"]
    ssibdry = nmtable["ssibdry"]
    isPlasma = nmtable["isPlasma"]
    btrz = nmtable["btrz"]
    brrz = nmtable["brrz"]
    bzrz = nmtable["bzrz"]
    bnorm = sqrt(btrz**2+brrz**2+bzrz**2)
    bp = sqrt(brrz**2+bzrz**2)
    bangle = arctan2(bp, bnorm)

    if viewer is None:
       from ifigure.interactive import figure
       plt = figure()
    else:
       plt =viewer
       plt.cls()

    plt.update(False)
       
    plt.subplot(3,3,[0,1,2], [3,4,5])
    plt.isec(3)
    plt.image(x,y,z, alpha=0.5)
    plt.contour(x,y,z,30)
    if len(nmtable["xlim"]) !=0:
       plt.plot(xlim, ylim)
    if len(nmtable["rbbbs"]) !=0:
       plt.plot(rb, zb)
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.title('PSI')

    plt.isec(4)
    plt.title('B-norm')
    plt.image(x,y, bnorm, alpha=0.5)
#    plt.image(x,y,bangle*180./3.1415, alpha=0.5)
    plt.contour(x,y,bnorm,30)
    if len(nmtable["xlim"]) !=0:
       plt.plot(xlim, ylim)
    if len(nmtable["rbbbs"]) !=0:
       plt.plot(rb, zb)
    plt.xlabel('R')
    plt.ylabel('Z')

    plt.isec(0)
    pres = nmtable["pres"]

    plt.plot(pres)
    plt.title('pressure')
    plt.xlabel('flux index')

    plt.isec(1)
    qpsi = nmtable["qpsi"]
    plt.plot(qpsi)
    plt.title('qpsi')
    plt.xlabel('flux index')

    plt.isec(2)
    fpol = nmtable["fpol"]
    plt.plot(fpol)
    plt.title('f_pol')
    plt.xlabel('flux index')
    plt.update(True)

    return plt


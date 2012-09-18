import math
from math import sqrt
from math import asin

import pychart
from pychart import *

import pylab
from pylab import *

import uuid

zero=-1
one=1

zero_h=0
one_h=1

setka=[i/10.0 for i in range(-10,10)]
setka_h=[i/10.0 for i in range(0,10)]
min_diff=10


h_l_glob=[zero_h,one_h]
l_glob=[zero,one]

results=[]


def polyFit2Var(xList,yList,ZList,degree=1):
	A = []
	n=len(xList)
	for i in range(n):
		A.append([])
		for xd in range(degree+1):
			for yd in range(degree+1-xd):
				A[i].append((xList[i]**xd)*(yList[i]**yd)) 

	c,_,_,_ = linalg.lstsq(A,ZList)
	j = 0
	for xd in range(0,degree+1):
		for yd in range(0,degree+1-xd):
			print " + (%.2f)*x^%d*y^%d" % (c[j], xd, yd),
			j += 1

def h_vrai(x):
    """
    if x>0:
        return float(1/(3.0-2*x))
    else:
        return float((1.0+x)/(3.0+x))
    """
    if x>0:
        return 1/(2-x)
    else:
        return (x+1)/(2+x)
    #return x
    #return math.sqrt(x)
    #return sin(x*math.pi/2.0)
    #return (exp(x)-1)/(exp(1)-1)
    
    
    
def h_1_vrai(x):
    #return x*x
    #return asin(x)*2.0/math.pi
    #return log(x*(exp(1)-1)+1)
    if x>0.5:
        return 2-1/x
    else:
        return (2*x-1)/(1-x)
    
def c(a,b):
   if a>=0 and b>=0:
       return a+b-a*b
   elif a<0 and b<0:
       return a+b+a*b
   else:
       return (a+b)/(1-min(abs(a),abs(b)))
   #return a*b

def d(a,b):
   return h_vrai(c(h_1_vrai(a),h_1_vrai(b)))
   #return (2*a*b)/((1-a)*(1-b)+2*a*b)
   #return a+b-a*b*(a+b-a*b)
   #return (a+b)/(1+a*b)
   #return 0.97-1.04*b  + 0.01*b*b  -0.53*a  + 2.07*a*b  - 0.56*a*a
  

def h(x):
    if x in l:
        return h_l[l.index(x)]
    else:
        l.append(x)
        h_l.append(x)
        return x

def is_monotonic(h_l):
    return h_l==sorted(h_l)

def polyFit(xList,yList,order=1):
    '''fit the data using a least squares and polynomial'''
    fList = [(lambda x,n=n: x**n) for n in range(order,-1,-1)]
    # build row for each element in y
    bList = []
    A_List = []
    for (thisX,thisY) in zip(xList,yList):
        bList.append(thisY)
        A_Row = [f(thisX) for f in fList]
        A_List.append(A_Row)
    b = matrix(bList).T
    A = matrix(A_List)
    w = inv(A.T*A)*A.T*b
    return w.T.tolist()[0]
    
def approx(w,x):
    if len(w)==3:
        return w[2]+w[1]*x+w[0]*x*x
    elif len(w)==4:
        return w[3] + w[2]*x+w[1]*x*x+w[0]*x*x*x
    elif len(w)==5:
        return w[4] + w[3]*x + w[2]*x*x+w[1]*x*x*x+w[0]*x*x*x*x
     
def is_monotonic(h_l,l,eps):
    is_good=True
    l_s=sorted(l)
    for i in range(1,len(l_s)):
        #it should be monotonic
        point_cur=l_s[i]
        index=l.index(point_cur)
        point_prev=l_s[i-1]
        index_prev=l.index(point_prev)
        if h_l[index] - h_l[index_prev] < -eps:
            is_good=False
            break
    return is_good

def is_not_plateau(h_l,l):
    is_good=True
    l_s=sorted(l)
    plateau=0
    for i in range(1,len(l_s)):
        #it should be monotonic
        point_cur=l_s[i]
        index=l.index(point_cur)
        point_prev=l_s[i-1]
        index_prev=l.index(point_prev)
        if abs(h_l[index] - h_l[index_prev]) ==0:
            plateau=plateau+1
          
    return float(plateau)/float(len(l) )    
 

def sequence(h_point_l,point_l, top,bottom, eps):
    l_c=[]
    h_l_c=[]
    h_point_l_cur=h_point_l
    point_l_cur=point_l
    l_c.append(point_l_cur)
    h_l_c.append(h_point_l_cur)
    if (point_l==0):
        return (l_c, h_l_c)
    
    while abs(point_l_cur-c(point_l_cur,point_l))>eps and abs( h_point_l_cur-d(h_point_l_cur,h_point_l))>eps:# check the condition
        point_l_cur=c(point_l_cur,point_l)
        h_point_l_cur=d(h_point_l_cur,h_point_l)
        l_c.append(point_l_cur)
        h_l_c.append(h_point_l_cur)
           
        
    if is_monotonic(h_l_c,l_c,0.1):
        return (l_c, h_l_c)
    else:
        return (0,0)
    


def add_sequence(h_l,l, top, bottom,eps, l_rest, min_diff,d_min):
    # l_rest - x points for which  there are no h(x) points defined
    # define them with the possible minimum diff and in accordance with the points already given
    #print len(results)
    if len(l_rest)==0:
        # READY
        #count diff
        #print "counting diff"
        w = polyFit(l,h_l,order=2)
        #if min(approx(w, -0.2)-approx(w, -0.3),approx(w, -0.1)-approx(w, -0.2),approx(w, 0.2)-approx(w, 0.1),approx(w, 0.3)-approx(w, 0.2))>0.06:
        max_diff=0
        for i in setka:
           try:
                max_diff=max_diff+abs(approx(w,c(i,i)) - d(approx(w,i),approx(w,i)))
           except ValueError:
                max_diff=max_diff+100# we don't need it
        
        max_diff=max_diff/100.0
        
        if  max_diff<1 and is_not_plateau(h_l,l)<1:
            min_diff=max_diff
            results.append(((l, h_l),max_diff))
            pylab.clf()
            pylab.xlabel(str(is_not_plateau(h_l,l)))
            pylab.ylabel(str(max_diff))
            #approx points
            pylab.plot(l,h_l, 'yo')
            #true points
            pylab.plot(sorted(l),[h_vrai(p) for p in sorted(l)], 'r-') 
            """            
            w = polyFit(l,[h_vrai(p) for p in l],order=3)
            #polynom for true points
            pylab.plot(sorted(l),[approx(w,p) for p in sorted(l)], 'g-') 
            """
            w = polyFit(l,h_l,order=3)
            #polynom for approx points
            pylab.plot(sorted(l),[approx(w,p) for p in sorted(l)], 'b-') 
            pylab.show()
            #pylab.savefig('./figs/Fig'+str(uuid.uuid4())+'.png')
            print max_diff
            #print 'f(x):='+str(w[0])+'*x^2+'+str(w[1])+'*x+'+str(w[2])+';'
            
            return None
           
            #print "added"
        return (l, h_l)
    else:
        point_l=l_rest[0]
        
        l_rest=l_rest[1:]
        
        bool_res=False
       
        for (h_point_l,diff) in d_min[point_l]:
            #print point_l, h_point_l 
            l_1=list(l)
            h_1=list(h_l)
            (seq_l,seq_h)=sequence(h_point_l,point_l, top,bottom, eps)
            if seq_l==0 and seq_h==0:
                continue
            h_1=h_1+seq_h
            l_1=l_1+seq_l  
            """                     
            h_setka=[]
            for p in setka[1:len(setka)]:
                try:
                    h_setka.append(h_1[l_1.index(p)])
                    
                except ValueError:
                    h_setka=h_setka
            """
            
            if  is_monotonic(h_1,l_1,0) and is_not_plateau(h_l,l)<0.3: 
                result=add_sequence(h_1,l_1, top, bottom,eps, l_rest, min_diff,d_min)
                if result==(l_1,h_1): # no changes, something went wrong
                    continue
                else:
                    bool_res=True
                    return result    
            else:
                """
                print "~~~~~~~~~~~~~~"
                if not is_monotonic(h_1,l_1,0):
                    print "because not monotonic", l_1,h_1

                if is_not_plateau(h_l,l)>=0.5:
                    print "because too many plateaus" 
                print "______________________"
                """
                continue
        # need to back track
        if not bool_res: 
            return (l, h_l)

            
def get_map():            
    d_min=dict([])
    for point_l in setka[1:len(setka)]:
        h=[]
        if point_l==0:
           d_min[point_l]=[(hpoint, 0) for hpoint in setka_h[1:len(setka_h)]] 
           continue
        for h_point_l in setka_h[1:len(setka_h)]:
            h_l=[zero_h,h_point_l,one_h]
            l=[zero,point_l,one]
            h_point_l_cur=h_point_l
            point_l_cur=point_l
            
            max_diff=0
           
            while one-point_l_cur>0.01 and point_l_cur-zero>0.01 and abs(point_l_cur-c(point_l_cur,point_l))>0.01:# check the condition
                
                point_l_cur=c(point_l_cur,point_l)
                h_point_l_cur=d(h_point_l_cur,h_point_l)
                l.append(point_l_cur)
                h_l.append(h_point_l_cur)
               
            w = polyFit(l,h_l,order=2)
                       
          
            for i in range(1, 100):
                point_i=i/100.0
                try:
                    max_diff=max_diff+abs(approx(w,c(point_i,point_i)) - d(approx(w,point_i),approx(w,point_i)))
                except ValueError:
                    max_diff=max_diff+100# we don't need it
            
            max_diff=max_diff/100.0
            i=0
            for (h_point,diff) in h:
                if diff>max_diff:
                    break
                i=i+1
            h.insert(i,(h_point_l,max_diff))
        d_min[point_l]=h 
    return d_min   
                


d_min=get_map()
"""
l=setka
pylab.plot(sorted(l),[h_vrai(p) for p in sorted(l)], 'r-') 
h_l=[zero_h]+[d_min[p][0][0] for p in setka[1:]]+[one_h]
w = polyFit(l,h_l,order=3)
#polynom for approx points
pylab.plot(sorted(l),[h_l[l.index(p)] for p in sorted(l)], 'bo') 
pylab.show()
print is_monotonic(h_l,l,0.01)
for p in l:
    print p, h_l[l.index(p)]
"""
add_sequence(h_l_glob,l_glob, one, zero, 0.1, setka[1:len(setka)],10,d_min)
"""


setka_1=[i/10.0 for i in range(0,10)]
setka_2=[]# 1111 then 2222 etc
setka_3=[]# 123 123 123 etc
for i in (1, len(setka_1)-2):
    setka_2=setka_2+[setka_1[i]]*len(setka_1) 
    setka_3=setka_3+setka_1
h_vrai_setka=[]
try:
   setka_4=[d(x,y) for (x,y) in zip(setka_2,setka_3)]
except ZeroDivisionError:
    print "s"

polyFit2Var(setka_2,setka_3,setka_4,2)

"""

        

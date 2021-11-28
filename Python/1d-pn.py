import re
import math
import os
import glob

def zerolistmaker(n):
    listofzeros = [0] * n
    return listofzeros

kk=zerolistmaker(398)
m=zerolistmaker(398)
N1=zerolistmaker(398)
N2=zerolistmaker(398)
K=zerolistmaker(398)
F=0

filename="3500K_0.55sifm-iprk-max.txt"
file=open(filename,"w")
file.close()

searchquery = '|psi|^2'
searchquery1 = '==== e('
searchquery2 = '     state #'
f=open('3500K_0.55sifm-pdos.out', 'r')
lines = f.readlines()                      #reading lines in f file
ar=[]
for j, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file
    if line.startswith(searchquery2):      #finding the lines containing searchquery2
       tr=lines[j]                         #naming lines j as tr
       er=re.findall("\d+", tr)            #extracting integer numbers in tr to a list of strings
       er1=map(int, er)                    #converting the list of strings to the list of numbers
       ar.append(er1)                      #putting the list to list ar
    if line.startswith(searchquery1):
       break
#print ar


for i, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file        
  a=[]
  a1=[]
  for M in range(1):                       # A trick for breaking if statement
    if line.startswith(searchquery1):      #finding the lines containing searchquery1
       t=lines[i]                          #naming lines i as t
       e=re.findall("\d+\.\d+", t)         #extracting decimal numbers in t to a list of strings
       e1=map(float, e)                    #converting the list of strings to the list of numbers
       e2=-sum(e1)                         #making all the numbers to negative numbers
       if "-" not in lines[i]:             #if not find "-" in the lines i   
          e2=sum(e1)                       #keep the numbers as positive numbers
       if not 5<e2<8:                   #only consider the state within this energy interval
          break
       F=F+1  
       for y in range(1, 100):            
           #if not line.startswith(searchquery):
           if "|psi|^2 =" not in lines[i+y]:         #if not find "|psi|^2 =" in lines i+y 
              s=lines[i+y]                           #naming lines i+y as s
              a.extend(re.findall("\d+\.\d+", s))    #extract all decimal numbers in line i+y 
              a1.extend(re.findall("\d+", s))        #putting integer numbers to a list of strings and append it to list "a1"
           else:
               b = map(float, a)                     #converting the list of string to the list of decimal numbers
               b1= map(int, a1)                      #converting the list of string to the list of integer numbers
               b2=[]
               for b2x in range(1, 1000):            #copying some numbers from list b1 to b2
                   k=3*b2x-1
                   if k<len(b1):
                      b2.append(b1[k])
                   else:
                      break
               #print b2
               #print b
               c=[x**4 for x in b]

               d=[x**2 for x in b]

               g=sum(d)*sum(d)
               
               for q in range(len(b2)):
                    N2[ar[b2[q]-1][1]-1]=N2[ar[b2[q]-1][1]-1]+(c[q]/g)              #the localization contribution of a given atom
                    K[ar[b2[q]-1][1]-1]=K[ar[b2[q]-1][1]-1]+1                     #the number of states located in a given atom
               for jjj in range(398):
                    if N1[jjj] < N2[jjj]:
                       N1[jjj] = N2[jjj]
               N2=zerolistmaker(398)
                    
               f=[x/g for x in c]

               g=sum(d)*sum(d)                           #calculating IPR value

               bb=sum(b)
               lb=len(b)
               ll=[]
               
               #print ll
               for ii in range(len(ll)):                   #printing the number of orbital, atom, orbital type, and localization
                    f1=open("3500K_0.55sifm-iprk-orbital.txt", "a")
                    f1.write(''+str(ar[b2[ii]-1][0])+' ')
                    f1.write(''+str(ar[b2[ii]-1][1])+' ')
                    f1.write(''+str(ar[b2[ii]-1][3])+' ')
                    f1.write(''+str(ar[b2[ii]-1][4])+' ')
                    f1.write(''+str(ll[ii])+'\n')
                    f1.close()
               for ii in range(len(ll)):
                   for jj in range(13):                           #dividing cell to small section of atom numberd
                      if 2*jj<ar[b2[ii]-1][1]<2*(jj+1.001):
                         #if ar[b2[ii]-1][3]==2:                   #calculating d orbitals in the localization
                            m[jj]=m[jj]+1                         # the number of d orbitals in creating localized states
                            #print m[jj]
                            kk[jj]=kk[jj]+ll[ii]                #the summation of localization

               #print e2, f
               #print f
               #print '*****'
               break
for jj in range(398):
   f1=open("3500K_0.55sifm-iprk-max.txt", "a")
   f1.write(''+str(jj+1)+' ')               #calculating z values.
   f1.write(''+str(N1[jj])+'\n')              #calculating average IPR
f1.close()
print (F)
#print f

import re
import math
import os
import glob

def zerolistmaker(n):
    listofzeros = [0] * n
    return listofzeros

kk=zerolistmaker(398)
m=zerolistmaker(398)
N1=zerolistmaker(398)
N2=zerolistmaker(398)
K=zerolistmaker(398)
F=0

filename="3500K_0.55sifm-iprk-max.txt"
file=open(filename,"w")
file.close()

searchquery = '|psi|^2'
searchquery1 = '==== e('
searchquery2 = '     state #'
f=open('3500K_0.55sifm-pdos.out', 'r')
lines = f.readlines()                      #reading lines in f file
ar=[]
for j, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file
    if line.startswith(searchquery2):      #finding the lines containing searchquery2
       tr=lines[j]                         #naming lines j as tr
       er=re.findall("\d+", tr)            #extracting integer numbers in tr to a list of strings
       er1=map(int, er)                    #converting the list of strings to the list of numbers
       ar.append(er1)                      #putting the list to list ar
    if line.startswith(searchquery1):
       break
#print ar


for i, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file        
  a=[]
  a1=[]
  for M in range(1):                       # A trick for breaking if statement
    if line.startswith(searchquery1):      #finding the lines containing searchquery1
       t=lines[i]                          #naming lines i as t
       e=re.findall("\d+\.\d+", t)         #extracting decimal numbers in t to a list of strings
       e1=map(float, e)                    #converting the list of strings to the list of numbers
       e2=-sum(e1)                         #making all the numbers to negative numbers
       if "-" not in lines[i]:             #if not find "-" in the lines i   
          e2=sum(e1)                       #keep the numbers as positive numbers
       if not 5<e2<8:                   #only consider the state within this energy interval
          break
       F=F+1  
       for y in range(1, 100):            
           #if not line.startswith(searchquery):
           if "|psi|^2 =" not in lines[i+y]:         #if not find "|psi|^2 =" in lines i+y 
              s=lines[i+y]                           #naming lines i+y as s
              a.extend(re.findall("\d+\.\d+", s))    #extract all decimal numbers in line i+y 
              a1.extend(re.findall("\d+", s))        #putting integer numbers to a list of strings and append it to list "a1"
           else:
               b = map(float, a)                     #converting the list of string to the list of decimal numbers
               b1= map(int, a1)                      #converting the list of string to the list of integer numbers
               b2=[]
               for b2x in range(1, 1000):            #copying some numbers from list b1 to b2
                   k=3*b2x-1
                   if k<len(b1):
                      b2.append(b1[k])
                   else:
                      break
               #print b2
               #print b
               c=[x**4 for x in b]

               d=[x**2 for x in b]

               g=sum(d)*sum(d)
               
               for q in range(len(b2)):
                    N2[ar[b2[q]-1][1]-1]=N2[ar[b2[q]-1][1]-1]+(c[q]/g)              #the localization contribution of a given atom
                    K[ar[b2[q]-1][1]-1]=K[ar[b2[q]-1][1]-1]+1                     #the number of states located in a given atom
               for jjj in range(398):
                    if N1[jjj] < N2[jjj]:
                       N1[jjj] = N2[jjj]
               N2=zerolistmaker(398)
                    
               f=[x/g for x in c]

               g=sum(d)*sum(d)                           #calculating IPR value

               bb=sum(b)
               lb=len(b)
               ll=[]
               
               #print ll
               for ii in range(len(ll)):                   #printing the number of orbital, atom, orbital type, and localization
                    f1=open("3500K_0.55sifm-iprk-orbital.txt", "a")
                    f1.write(''+str(ar[b2[ii]-1][0])+' ')
                    f1.write(''+str(ar[b2[ii]-1][1])+' ')
                    f1.write(''+str(ar[b2[ii]-1][3])+' ')
                    f1.write(''+str(ar[b2[ii]-1][4])+' ')
                    f1.write(''+str(ll[ii])+'\n')
                    f1.close()
               for ii in range(len(ll)):
                   for jj in range(13):                           #dividing cell to small section of atom numberd
                      if 2*jj<ar[b2[ii]-1][1]<2*(jj+1.001):
                         #if ar[b2[ii]-1][3]==2:                   #calculating d orbitals in the localization
                            m[jj]=m[jj]+1                         # the number of d orbitals in creating localized states
                            #print m[jj]
                            kk[jj]=kk[jj]+ll[ii]                #the summation of localization

               #print e2, f
               #print f
               #print '*****'
               break
for jj in range(398):
   f1=open("3500K_0.55sifm-iprk-max.txt", "a")
   f1.write(''+str(jj+1)+' ')               #calculating z values.
   f1.write(''+str(N1[jj])+'\n')              #calculating average IPR
f1.close()
print (F)
#print f


               

import re
import math
import os
import glob

def zerolistmaker(n):
    listofzeros = [0] * n
    return listofzeros

kk=zerolistmaker(398)
m=zerolistmaker(398)
N1=zerolistmaker(398)
N2=zerolistmaker(398)
K=zerolistmaker(398)
F=0

filename="3500K_0.55sifm-iprk-max.txt"
file=open(filename,"w")
file.close()

searchquery = '|psi|^2'
searchquery1 = '==== e('
searchquery2 = '     state #'
f=open('3500K_0.55sifm-pdos.out', 'r')
lines = f.readlines()                      #reading lines in f file
ar=[]
for j, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file
    if line.startswith(searchquery2):      #finding the lines containing searchquery2
       tr=lines[j]                         #naming lines j as tr
       er=re.findall("\d+", tr)            #extracting integer numbers in tr to a list of strings
       er1=map(int, er)                    #converting the list of strings to the list of numbers
       ar.append(er1)                      #putting the list to list ar
    if line.startswith(searchquery1):
       break
#print ar


for i, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file        
  a=[]
  a1=[]
  for M in range(1):                       # A trick for breaking if statement
    if line.startswith(searchquery1):      #finding the lines containing searchquery1
       t=lines[i]                          #naming lines i as t
       e=re.findall("\d+\.\d+", t)         #extracting decimal numbers in t to a list of strings
       e1=map(float, e)                    #converting the list of strings to the list of numbers
       e2=-sum(e1)                         #making all the numbers to negative numbers
       if "-" not in lines[i]:             #if not find "-" in the lines i   
          e2=sum(e1)                       #keep the numbers as positive numbers
       if not 5<e2<8:                   #only consider the state within this energy interval
          break
       F=F+1  
       for y in range(1, 100):            
           #if not line.startswith(searchquery):
           if "|psi|^2 =" not in lines[i+y]:         #if not find "|psi|^2 =" in lines i+y 
              s=lines[i+y]                           #naming lines i+y as s
              a.extend(re.findall("\d+\.\d+", s))    #extract all decimal numbers in line i+y 
              a1.extend(re.findall("\d+", s))        #putting integer numbers to a list of strings and append it to list "a1"
           else:
               b = map(float, a)                     #converting the list of string to the list of decimal numbers
               b1= map(int, a1)                      #converting the list of string to the list of integer numbers
               b2=[]
               for b2x in range(1, 1000):            #copying some numbers from list b1 to b2
                   k=3*b2x-1
                   if k<len(b1):
                      b2.append(b1[k])
                   else:
                      break
               #print b2
               #print b
               c=[x**4 for x in b]

               d=[x**2 for x in b]

               g=sum(d)*sum(d)
               
               for q in range(len(b2)):
                    N2[ar[b2[q]-1][1]-1]=N2[ar[b2[q]-1][1]-1]+(c[q]/g)              #the localization contribution of a given atom
                    K[ar[b2[q]-1][1]-1]=K[ar[b2[q]-1][1]-1]+1                     #the number of states located in a given atom
               for jjj in range(398):
                    if N1[jjj] < N2[jjj]:
                       N1[jjj] = N2[jjj]
               N2=zerolistmaker(398)
                    
               f=[x/g for x in c]

               g=sum(d)*sum(d)                           #calculating IPR value

               bb=sum(b)
               lb=len(b)
               ll=[]
               
               #print ll
               for ii in range(len(ll)):                   #printing the number of orbital, atom, orbital type, and localization
                    f1=open("3500K_0.55sifm-iprk-orbital.txt", "a")
                    f1.write(''+str(ar[b2[ii]-1][0])+' ')
                    f1.write(''+str(ar[b2[ii]-1][1])+' ')
                    f1.write(''+str(ar[b2[ii]-1][3])+' ')
                    f1.write(''+str(ar[b2[ii]-1][4])+' ')
                    f1.write(''+str(ll[ii])+'\n')
                    f1.close()
               for ii in range(len(ll)):
                   for jj in range(13):                           #dividing cell to small section of atom numberd
                      if 2*jj<ar[b2[ii]-1][1]<2*(jj+1.001):
                         #if ar[b2[ii]-1][3]==2:                   #calculating d orbitals in the localization
                            m[jj]=m[jj]+1                         # the number of d orbitals in creating localized states
                            #print m[jj]
                            kk[jj]=kk[jj]+ll[ii]                #the summation of localization

               #print e2, f
               #print f
               #print '*****'
               break
for jj in range(398):
   f1=open("3500K_0.55sifm-iprk-max.txt", "a")
   f1.write(''+str(jj+1)+' ')               #calculating z values.
   f1.write(''+str(N1[jj])+'\n')              #calculating average IPR
f1.close()
print (F)
#print f


               

import re
import math
import os
import glob

def zerolistmaker(n):
    listofzeros = [0] * n
    return listofzeros

kk=zerolistmaker(398)
m=zerolistmaker(398)
N1=zerolistmaker(398)
N2=zerolistmaker(398)
K=zerolistmaker(398)
F=0

filename="3500K_0.55sifm-iprk-max.txt"
file=open(filename,"w")
file.close()

searchquery = '|psi|^2'
searchquery1 = '==== e('
searchquery2 = '     state #'
f=open('3500K_0.55sifm-pdos.out', 'r')
lines = f.readlines()                      #reading lines in f file
ar=[]
for j, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file
    if line.startswith(searchquery2):      #finding the lines containing searchquery2
       tr=lines[j]                         #naming lines j as tr
       er=re.findall("\d+", tr)            #extracting integer numbers in tr to a list of strings
       er1=map(int, er)                    #converting the list of strings to the list of numbers
       ar.append(er1)                      #putting the list to list ar
    if line.startswith(searchquery1):
       break
#print ar


for i, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file        
  a=[]
  a1=[]
  for M in range(1):                       # A trick for breaking if statement
    if line.startswith(searchquery1):      #finding the lines containing searchquery1
       t=lines[i]                          #naming lines i as t
       e=re.findall("\d+\.\d+", t)         #extracting decimal numbers in t to a list of strings
       e1=map(float, e)                    #converting the list of strings to the list of numbers
       e2=-sum(e1)                         #making all the numbers to negative numbers
       if "-" not in lines[i]:             #if not find "-" in the lines i   
          e2=sum(e1)                       #keep the numbers as positive numbers
       if not 5<e2<8:                   #only consider the state within this energy interval
          break
       F=F+1  
       for y in range(1, 100):            
           #if not line.startswith(searchquery):
           if "|psi|^2 =" not in lines[i+y]:         #if not find "|psi|^2 =" in lines i+y 
              s=lines[i+y]                           #naming lines i+y as s
              a.extend(re.findall("\d+\.\d+", s))    #extract all decimal numbers in line i+y 
              a1.extend(re.findall("\d+", s))        #putting integer numbers to a list of strings and append it to list "a1"
           else:
               b = map(float, a)                     #converting the list of string to the list of decimal numbers
               b1= map(int, a1)                      #converting the list of string to the list of integer numbers
               b2=[]
               for b2x in range(1, 1000):            #copying some numbers from list b1 to b2
                   k=3*b2x-1
                   if k<len(b1):
                      b2.append(b1[k])
                   else:
                      break
               #print b2
               #print b
               c=[x**4 for x in b]

               d=[x**2 for x in b]

               g=sum(d)*sum(d)
               
               for q in range(len(b2)):
                    N2[ar[b2[q]-1][1]-1]=N2[ar[b2[q]-1][1]-1]+(c[q]/g)              #the localization contribution of a given atom
                    K[ar[b2[q]-1][1]-1]=K[ar[b2[q]-1][1]-1]+1                     #the number of states located in a given atom
               for jjj in range(398):
                    if N1[jjj] < N2[jjj]:
                       N1[jjj] = N2[jjj]
               N2=zerolistmaker(398)
                    
               f=[x/g for x in c]

               g=sum(d)*sum(d)                           #calculating IPR value

               bb=sum(b)
               lb=len(b)
               ll=[]
               
               #print ll
               for ii in range(len(ll)):                   #printing the number of orbital, atom, orbital type, and localization
                    f1=open("3500K_0.55sifm-iprk-orbital.txt", "a")
                    f1.write(''+str(ar[b2[ii]-1][0])+' ')
                    f1.write(''+str(ar[b2[ii]-1][1])+' ')
                    f1.write(''+str(ar[b2[ii]-1][3])+' ')
                    f1.write(''+str(ar[b2[ii]-1][4])+' ')
                    f1.write(''+str(ll[ii])+'\n')
                    f1.close()
               for ii in range(len(ll)):
                   for jj in range(13):                           #dividing cell to small section of atom numberd
                      if 2*jj<ar[b2[ii]-1][1]<2*(jj+1.001):
                         #if ar[b2[ii]-1][3]==2:                   #calculating d orbitals in the localization
                            m[jj]=m[jj]+1                         # the number of d orbitals in creating localized states
                            #print m[jj]
                            kk[jj]=kk[jj]+ll[ii]                #the summation of localization

               #print e2, f
               #print f
               #print '*****'
               break
for jj in range(398):
   f1=open("3500K_0.55sifm-iprk-max.txt", "a")
   f1.write(''+str(jj+1)+' ')               #calculating z values.
   f1.write(''+str(N1[jj])+'\n')              #calculating average IPR
f1.close()
print (F)
#print f


               

import re
import math
import os
import glob

def zerolistmaker(n):
    listofzeros = [0] * n
    return listofzeros

kk=zerolistmaker(398)
m=zerolistmaker(398)
N1=zerolistmaker(398)
N2=zerolistmaker(398)
K=zerolistmaker(398)
F=0

filename="3500K_0.55sifm-iprk-max.txt"
file=open(filename,"w")
file.close()

searchquery = '|psi|^2'
searchquery1 = '==== e('
searchquery2 = '     state #'
f=open('3500K_0.55sifm-pdos.out', 'r')
lines = f.readlines()                      #reading lines in f file
ar=[]
for j, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file
    if line.startswith(searchquery2):      #finding the lines containing searchquery2
       tr=lines[j]                         #naming lines j as tr
       er=re.findall("\d+", tr)            #extracting integer numbers in tr to a list of strings
       er1=map(int, er)                    #converting the list of strings to the list of numbers
       ar.append(er1)                      #putting the list to list ar
    if line.startswith(searchquery1):
       break
#print ar


for i, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file        
  a=[]
  a1=[]
  for M in range(1):                       # A trick for breaking if statement
    if line.startswith(searchquery1):      #finding the lines containing searchquery1
       t=lines[i]                          #naming lines i as t
       e=re.findall("\d+\.\d+", t)         #extracting decimal numbers in t to a list of strings
       e1=map(float, e)                    #converting the list of strings to the list of numbers
       e2=-sum(e1)                         #making all the numbers to negative numbers
       if "-" not in lines[i]:             #if not find "-" in the lines i   
          e2=sum(e1)                       #keep the numbers as positive numbers
       if not 5<e2<8:                   #only consider the state within this energy interval
          break
       F=F+1  
       for y in range(1, 100):            
           #if not line.startswith(searchquery):
           if "|psi|^2 =" not in lines[i+y]:         #if not find "|psi|^2 =" in lines i+y 
              s=lines[i+y]                           #naming lines i+y as s
              a.extend(re.findall("\d+\.\d+", s))    #extract all decimal numbers in line i+y 
              a1.extend(re.findall("\d+", s))        #putting integer numbers to a list of strings and append it to list "a1"
           else:
               b = map(float, a)                     #converting the list of string to the list of decimal numbers
               b1= map(int, a1)                      #converting the list of string to the list of integer numbers
               b2=[]
               for b2x in range(1, 1000):            #copying some numbers from list b1 to b2
                   k=3*b2x-1
                   if k<len(b1):
                      b2.append(b1[k])
                   else:
                      break
               #print b2
               #print b
               c=[x**4 for x in b]

               d=[x**2 for x in b]

               g=sum(d)*sum(d)
               
               for q in range(len(b2)):
                    N2[ar[b2[q]-1][1]-1]=N2[ar[b2[q]-1][1]-1]+(c[q]/g)              #the localization contribution of a given atom
                    K[ar[b2[q]-1][1]-1]=K[ar[b2[q]-1][1]-1]+1                     #the number of states located in a given atom
               for jjj in range(398):
                    if N1[jjj] < N2[jjj]:
                       N1[jjj] = N2[jjj]
               N2=zerolistmaker(398)
                    
               f=[x/g for x in c]

               g=sum(d)*sum(d)                           #calculating IPR value

               bb=sum(b)
               lb=len(b)
               ll=[]
               
               #print ll
               for ii in range(len(ll)):                   #printing the number of orbital, atom, orbital type, and localization
                    f1=open("3500K_0.55sifm-iprk-orbital.txt", "a")
                    f1.write(''+str(ar[b2[ii]-1][0])+' ')
                    f1.write(''+str(ar[b2[ii]-1][1])+' ')
                    f1.write(''+str(ar[b2[ii]-1][3])+' ')
                    f1.write(''+str(ar[b2[ii]-1][4])+' ')
                    f1.write(''+str(ll[ii])+'\n')
                    f1.close()
               for ii in range(len(ll)):
                   for jj in range(13):                           #dividing cell to small section of atom numberd
                      if 2*jj<ar[b2[ii]-1][1]<2*(jj+1.001):
                         #if ar[b2[ii]-1][3]==2:                   #calculating d orbitals in the localization
                            m[jj]=m[jj]+1                         # the number of d orbitals in creating localized states
                            #print m[jj]
                            kk[jj]=kk[jj]+ll[ii]                #the summation of localization

               #print e2, f
               #print f
               #print '*****'
               break
for jj in range(398):
   f1=open("3500K_0.55sifm-iprk-max.txt", "a")
   f1.write(''+str(jj+1)+' ')               #calculating z values.
   f1.write(''+str(N1[jj])+'\n')              #calculating average IPR
f1.close()
print (F)
#print f


               

import re
import math
import os
import glob

def zerolistmaker(n):
    listofzeros = [0] * n
    return listofzeros

kk=zerolistmaker(398)
m=zerolistmaker(398)
N1=zerolistmaker(398)
N2=zerolistmaker(398)
K=zerolistmaker(398)
F=0

filename="3500K_0.55sifm-iprk-max.txt"
file=open(filename,"w")
file.close()

searchquery = '|psi|^2'
searchquery1 = '==== e('
searchquery2 = '     state #'
f=open('3500K_0.55sifm-pdos.out', 'r')
lines = f.readlines()                      #reading lines in f file
ar=[]
for j, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file
    if line.startswith(searchquery2):      #finding the lines containing searchquery2
       tr=lines[j]                         #naming lines j as tr
       er=re.findall("\d+", tr)            #extracting integer numbers in tr to a list of strings
       er1=map(int, er)                    #converting the list of strings to the list of numbers
       ar.append(er1)                      #putting the list to list ar
    if line.startswith(searchquery1):
       break
#print ar


for i, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file        
  a=[]
  a1=[]
  for M in range(1):                       # A trick for breaking if statement
    if line.startswith(searchquery1):      #finding the lines containing searchquery1
       t=lines[i]                          #naming lines i as t
       e=re.findall("\d+\.\d+", t)         #extracting decimal numbers in t to a list of strings
       e1=map(float, e)                    #converting the list of strings to the list of numbers
       e2=-sum(e1)                         #making all the numbers to negative numbers
       if "-" not in lines[i]:             #if not find "-" in the lines i   
          e2=sum(e1)                       #keep the numbers as positive numbers
       if not 5<e2<8:                   #only consider the state within this energy interval
          break
       F=F+1  
       for y in range(1, 100):            
           #if not line.startswith(searchquery):
           if "|psi|^2 =" not in lines[i+y]:         #if not find "|psi|^2 =" in lines i+y 
              s=lines[i+y]                           #naming lines i+y as s
              a.extend(re.findall("\d+\.\d+", s))    #extract all decimal numbers in line i+y 
              a1.extend(re.findall("\d+", s))        #putting integer numbers to a list of strings and append it to list "a1"
           else:
               b = map(float, a)                     #converting the list of string to the list of decimal numbers
               b1= map(int, a1)                      #converting the list of string to the list of integer numbers
               b2=[]
               for b2x in range(1, 1000):            #copying some numbers from list b1 to b2
                   k=3*b2x-1
                   if k<len(b1):
                      b2.append(b1[k])
                   else:
                      break
               #print b2
               #print b
               c=[x**4 for x in b]

               d=[x**2 for x in b]

               g=sum(d)*sum(d)
               
               for q in range(len(b2)):
                    N2[ar[b2[q]-1][1]-1]=N2[ar[b2[q]-1][1]-1]+(c[q]/g)              #the localization contribution of a given atom
                    K[ar[b2[q]-1][1]-1]=K[ar[b2[q]-1][1]-1]+1                     #the number of states located in a given atom
               for jjj in range(398):
                    if N1[jjj] < N2[jjj]:
                       N1[jjj] = N2[jjj]
               N2=zerolistmaker(398)
                    
               f=[x/g for x in c]

               g=sum(d)*sum(d)                           #calculating IPR value

               bb=sum(b)
               lb=len(b)
               ll=[]
               
               #print ll
               for ii in range(len(ll)):                   #printing the number of orbital, atom, orbital type, and localization
                    f1=open("3500K_0.55sifm-iprk-orbital.txt", "a")
                    f1.write(''+str(ar[b2[ii]-1][0])+' ')
                    f1.write(''+str(ar[b2[ii]-1][1])+' ')
                    f1.write(''+str(ar[b2[ii]-1][3])+' ')
                    f1.write(''+str(ar[b2[ii]-1][4])+' ')
                    f1.write(''+str(ll[ii])+'\n')
                    f1.close()
               for ii in range(len(ll)):
                   for jj in range(13):                           #dividing cell to small section of atom numberd
                      if 2*jj<ar[b2[ii]-1][1]<2*(jj+1.001):
                         #if ar[b2[ii]-1][3]==2:                   #calculating d orbitals in the localization
                            m[jj]=m[jj]+1                         # the number of d orbitals in creating localized states
                            #print m[jj]
                            kk[jj]=kk[jj]+ll[ii]                #the summation of localization

               #print e2, f
               #print f
               #print '*****'
               break
for jj in range(398):
   f1=open("3500K_0.55sifm-iprk-max.txt", "a")
   f1.write(''+str(jj+1)+' ')               #calculating z values.
   f1.write(''+str(N1[jj])+'\n')              #calculating average IPR
f1.close()
print (F)
#print f


               

import re
import math
import os
import glob

def zerolistmaker(n):
    listofzeros = [0] * n
    return listofzeros

kk=zerolistmaker(398)
m=zerolistmaker(398)
N1=zerolistmaker(398)
N2=zerolistmaker(398)
K=zerolistmaker(398)
F=0

filename="3500K_0.55sifm-iprk-max.txt"
file=open(filename,"w")
file.close()

searchquery = '|psi|^2'
searchquery1 = '==== e('
searchquery2 = '     state #'
f=open('3500K_0.55sifm-pdos.out', 'r')
lines = f.readlines()                      #reading lines in f file
ar=[]
for j, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file
    if line.startswith(searchquery2):      #finding the lines containing searchquery2
       tr=lines[j]                         #naming lines j as tr
       er=re.findall("\d+", tr)            #extracting integer numbers in tr to a list of strings
       er1=map(int, er)                    #converting the list of strings to the list of numbers
       ar.append(er1)                      #putting the list to list ar
    if line.startswith(searchquery1):
       break
#print ar


for i, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file        
  a=[]
  a1=[]
  for M in range(1):                       # A trick for breaking if statement
    if line.startswith(searchquery1):      #finding the lines containing searchquery1
       t=lines[i]                          #naming lines i as t
       e=re.findall("\d+\.\d+", t)         #extracting decimal numbers in t to a list of strings
       e1=map(float, e)                    #converting the list of strings to the list of numbers
       e2=-sum(e1)                         #making all the numbers to negative numbers
       if "-" not in lines[i]:             #if not find "-" in the lines i   
          e2=sum(e1)                       #keep the numbers as positive numbers
       if not 5<e2<8:                   #only consider the state within this energy interval
          break
       F=F+1  
       for y in range(1, 100):            
           #if not line.startswith(searchquery):
           if "|psi|^2 =" not in lines[i+y]:         #if not find "|psi|^2 =" in lines i+y 
              s=lines[i+y]                           #naming lines i+y as s
              a.extend(re.findall("\d+\.\d+", s))    #extract all decimal numbers in line i+y 
              a1.extend(re.findall("\d+", s))        #putting integer numbers to a list of strings and append it to list "a1"
           else:
               b = map(float, a)                     #converting the list of string to the list of decimal numbers
               b1= map(int, a1)                      #converting the list of string to the list of integer numbers
               b2=[]
               for b2x in range(1, 1000):            #copying some numbers from list b1 to b2
                   k=3*b2x-1
                   if k<len(b1):
                      b2.append(b1[k])
                   else:
                      break
               #print b2
               #print b
               c=[x**4 for x in b]

               d=[x**2 for x in b]

               g=sum(d)*sum(d)
               
               for q in range(len(b2)):
                    N2[ar[b2[q]-1][1]-1]=N2[ar[b2[q]-1][1]-1]+(c[q]/g)              #the localization contribution of a given atom
                    K[ar[b2[q]-1][1]-1]=K[ar[b2[q]-1][1]-1]+1                     #the number of states located in a given atom
               for jjj in range(398):
                    if N1[jjj] < N2[jjj]:
                       N1[jjj] = N2[jjj]
               N2=zerolistmaker(398)
                    
               f=[x/g for x in c]

               g=sum(d)*sum(d)                           #calculating IPR value

               bb=sum(b)
               lb=len(b)
               ll=[]
               
               #print ll
               for ii in range(len(ll)):                   #printing the number of orbital, atom, orbital type, and localization
                    f1=open("3500K_0.55sifm-iprk-orbital.txt", "a")
                    f1.write(''+str(ar[b2[ii]-1][0])+' ')
                    f1.write(''+str(ar[b2[ii]-1][1])+' ')
                    f1.write(''+str(ar[b2[ii]-1][3])+' ')
                    f1.write(''+str(ar[b2[ii]-1][4])+' ')
                    f1.write(''+str(ll[ii])+'\n')
                    f1.close()
               for ii in range(len(ll)):
                   for jj in range(13):                           #dividing cell to small section of atom numberd
                      if 2*jj<ar[b2[ii]-1][1]<2*(jj+1.001):
                         #if ar[b2[ii]-1][3]==2:                   #calculating d orbitals in the localization
                            m[jj]=m[jj]+1                         # the number of d orbitals in creating localized states
                            #print m[jj]
                            kk[jj]=kk[jj]+ll[ii]                #the summation of localization

               #print e2, f
               #print f
               #print '*****'
               break
for jj in range(398):
   f1=open("3500K_0.55sifm-iprk-max.txt", "a")
   f1.write(''+str(jj+1)+' ')               #calculating z values.
   f1.write(''+str(N1[jj])+'\n')              #calculating average IPR
f1.close()
print (F)
#print f


               

import re
import math
import os
import glob

def zerolistmaker(n):
    listofzeros = [0] * n
    return listofzeros

kk=zerolistmaker(398)
m=zerolistmaker(398)
N1=zerolistmaker(398)
N2=zerolistmaker(398)
K=zerolistmaker(398)
F=0

filename="3500K_0.55sifm-iprk-max.txt"
file=open(filename,"w")
file.close()

searchquery = '|psi|^2'
searchquery1 = '==== e('
searchquery2 = '     state #'
f=open('3500K_0.55sifm-pdos.out', 'r')
lines = f.readlines()                      #reading lines in f file
ar=[]
for j, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file
    if line.startswith(searchquery2):      #finding the lines containing searchquery2
       tr=lines[j]                         #naming lines j as tr
       er=re.findall("\d+", tr)            #extracting integer numbers in tr to a list of strings
       er1=map(int, er)                    #converting the list of strings to the list of numbers
       ar.append(er1)                      #putting the list to list ar
    if line.startswith(searchquery1):
       break
#print ar


for i, line in enumerate(lines):           #counting lines and do the following jobs for all lines in the file        
  a=[]
  a1=[]
  for M in range(1):                       # A trick for breaking if statement
    if line.startswith(searchquery1):      #finding the lines containing searchquery1
       t=lines[i]                          #naming lines i as t
       e=re.findall("\d+\.\d+", t)         #extracting decimal numbers in t to a list of strings
       e1=map(float, e)                    #converting the list of strings to the list of numbers
       e2=-sum(e1)                         #making all the numbers to negative numbers
       if "-" not in lines[i]:             #if not find "-" in the lines i   
          e2=sum(e1)                       #keep the numbers as positive numbers
       if not 5<e2<8:                   #only consider the state within this energy interval
          break
       F=F+1  
       for y in range(1, 100):            
           #if not line.startswith(searchquery):
           if "|psi|^2 =" not in lines[i+y]:         #if not find "|psi|^2 =" in lines i+y 
              s=lines[i+y]                           #naming lines i+y as s
              a.extend(re.findall("\d+\.\d+", s))    #extract all decimal numbers in line i+y 
              a1.extend(re.findall("\d+", s))        #putting integer numbers to a list of strings and append it to list "a1"
           else:
               b = map(float, a)                     #converting the list of string to the list of decimal numbers
               b1= map(int, a1)                      #converting the list of string to the list of integer numbers
               b2=[]
               for b2x in range(1, 1000):            #copying some numbers from list b1 to b2
                   k=3*b2x-1
                   if k<len(b1):
                      b2.append(b1[k])
                   else:
                      break
               #print b2
               #print b
               c=[x**4 for x in b]

               d=[x**2 for x in b]

               g=sum(d)*sum(d)
               
               for q in range(len(b2)):
                    N2[ar[b2[q]-1][1]-1]=N2[ar[b2[q]-1][1]-1]+(c[q]/g)              #the localization contribution of a given atom
                    K[ar[b2[q]-1][1]-1]=K[ar[b2[q]-1][1]-1]+1                     #the number of states located in a given atom
               for jjj in range(398):
                    if N1[jjj] < N2[jjj]:
                       N1[jjj] = N2[jjj]
               N2=zerolistmaker(398)
                    
               f=[x/g for x in c]

               g=sum(d)*sum(d)                           #calculating IPR value

               bb=sum(b)
               lb=len(b)
               ll=[]
               
               #print ll
               for ii in range(len(ll)):                   #printing the number of orbital, atom, orbital type, and localization
                    f1=open("3500K_0.55sifm-iprk-orbital.txt", "a")
                    f1.write(''+str(ar[b2[ii]-1][0])+' ')
                    f1.write(''+str(ar[b2[ii]-1][1])+' ')
                    f1.write(''+str(ar[b2[ii]-1][3])+' ')
                    f1.write(''+str(ar[b2[ii]-1][4])+' ')
                    f1.write(''+str(ll[ii])+'\n')
                    f1.close()
               for ii in range(len(ll)):
                   for jj in range(13):                           #dividing cell to small section of atom numberd
                      if 2*jj<ar[b2[ii]-1][1]<2*(jj+1.001):
                         #if ar[b2[ii]-1][3]==2:                   #calculating d orbitals in the localization
                            m[jj]=m[jj]+1                         # the number of d orbitals in creating localized states
                            #print m[jj]
                            kk[jj]=kk[jj]+ll[ii]                #the summation of localization

               #print e2, f
               #print f
               #print '*****'
               break
for jj in range(398):
   f1=open("3500K_0.55sifm-iprk-max.txt", "a")
   f1.write(''+str(jj+1)+' ')               #calculating z values.
   f1.write(''+str(N1[jj])+'\n')              #calculating average IPR
f1.close()
print (F)
#print f


               




               



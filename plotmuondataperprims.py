import numpy as np
from scipy.special import kn
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import math
import matplotlib.ticker
from matplotlib.ticker import (MultipleLocator, 
                               FormatStrFormatter, 
                               AutoMinorLocator) 
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
k = 1.38064852*10**(-23)
c=3*10**8
m=1.883531627*10**(-28)
path="/Users/lukecalvin/2021/chicane_sims/outof2cm/finished_files copy/"
#i = input("-- specify particle type --\nenter number-\nmuons: 26\npositrons: 27\nelecrons: 28\nphotons: 29\n:" )
#search_phrase = input("-- specify data set --" )
i=27
search_phrase ="gamma_Ol"
#q = open('muon_production_%d_sum.lis'%(int(i)),'r')
#line1 = q.readlines()[4]
#parts = line1.split(' ')    
prim = 1000000000
#print(prim)
#prim = input("enter number of primaries:")
total = 0
#logx = input("log plot x axis?\n0:no\n1:yes\n:")
#logy = input("log plot y axis?\n0:no\n1:yes\n:")
#prim = 93000000
error=[]
errorp=[]
errore=[]
errorg=[]
f = open('%smuon_production_%d_tab.lis'%(path,int(i)),'r')
x=[]
prob=[]
gamma=[]
K2=[]
p=[]
xs=[]
y=[]
yg=[]
xg=[]
ys=[]
errors=[]
ymax=0
norm=0.1
normmj=720962.4432
maxprob=0
print(xs)
line_num = 0
line_num_2 = 0
end_line = 100000000000000000000000000
start_line = 1000000000000
for line in f.readlines():
    index=0
    end_count = 0
    if line.find(search_phrase) >= 0:
        start_line=line_num+1
       # print(line_num)
        #print(line)
   # print("line_num:",line_num)
   # print("start_line:",start_line)
    if line_num > start_line:
         if line.find(" ", 2,3) >= 0:  
            end_line=line_num
            #print("end line:", end_line)
    if line_num > start_line and line_num < end_line:
        #print(line)
        check = (line.split())
        if len(check) == 0:
            break
        ys = float(line.split()[2])
        errors = float(line.split()[3])
        if ys != 0:
            if ymax < ys*norm:
                ymax = ys*norm
            y.append(ys*norm)
            yg.append(ys*norm)
            errorg.append((ys*norm/100)*errors)
            total=total+(ys*norm)
            xs = float(line.split()[0])
            
            x.append(xs)
            xg.append(xs)
            #print("xs: ",xs)
            if ((((1.6*10**(-19))*1000000000*xs)**2)-((m*(c**2))**2))>0:
                p.append((math.sqrt((((1.6*10**(-19))*1000000000*xs)**2)-((m*(c**2))**2)))/c)
            if ((((1.6*10**(-19))*1000000000*xs)**2)-((m*(c**2))**2))<0:
                p.append(0)
            #print("p: ",p)
            K2.append(0)
            prob.append(0)
            gamma.append(0)
        #print(x)
    if line_num >= end_line:
        break
    line_num += 1
totalprob=0
k = 1.38064852*10**(-23)
c=3*10**8
m=1.883531627*10**(-28)
lowchisq=1000000000000000000000000000000
#T=[]
#for s in range(0,600,50):
#    T.append(s)
testprob=[]
testy=[]
testerror=[]
for z in range(len(p)):
    testprob.append(0)
    testy.append(0)
    testerror.append(0)
for T in range(8000000000000,8830000000000,1000000000):
    testprob=[]
    testy=[]
    testerror=[]
    for z in range(len(p)):
        testprob.append(0)
        testy.append(0)
        testerror.append(0)
    totalprob=0
    maxprob=0
    print("TEMPERATURE",T)
    theta = k*T/(m*(c**2))
    #print("theta: ",theta)
    for z in range(len(p)):
        gamma[z]=math.sqrt(1+(p[z]/(m*c))**2)
        #print("gamma: ",gamma)
        K2[z] = kn(2,1/theta)
        #K2[z]=((p[z]-1)/((gamma[0]**(-p[z]+1))-(gamma[len(gamma)-1]**(-p[z]+1))))*(gamma[z]**(-p[z]))
        #print("K2: ",K2)
        prob[z]=((((gamma[z]**2)*(math.sqrt(1-(1/(gamma[z]**2)))))/(theta*K2[z]))*math.exp(-gamma[z]/theta))
        totalprob=totalprob+prob[z]
        if prob[z]>maxprob:
            maxprob=prob[z]
        #print("prob: ",prob)
    for z in range(len(p)):
        testprob[z]=prob[z]/totalprob
        testy[z]=y[z]/total
        testerror[z]=errorg[z]/total
    summ=0
    v=len(p)-1-2
    for z in range(len(p)):
        summ=summ+(((testy[z]-testprob[z])/testerror[z])**2)
    chisq=(1/v)*summ
    print("chisq: ",chisq)
    if chisq<lowchisq:
        lowchisq=chisq
        optT=T
        print("optimal T: ",optT)
T=optT
totalprob=0
maxprob=0
theta = k*T/(m*(c**2))
print("theta: ",theta)
for z in range(len(p)):
    gamma[z]=math.sqrt(1+(p[z]/(m*c))**2)
    print("gamma: ",gamma)
    K2[z] = kn(2,1/theta)
    #K2[z]=((p[z]-1)/((gamma[0]**(-p[z]+1))-(gamma[len(gamma)-1]**(-p[z]+1))))*(gamma[z]**(-p[z]))
    print("K2: ",K2)
    prob[z]=((((gamma[z]**2)*(math.sqrt(1-(1/(gamma[z]**2)))))/(theta*K2[z]))*math.exp(-gamma[z]/theta))
    totalprob=totalprob+prob[z]
    if prob[z]>maxprob:
        maxprob=prob[z]
    print("prob: ",prob)
for z in range(len(p)):
    prob[z]=(prob[z]/totalprob)*total
    y[z]=y[z]
#    error[z]=error[z]
summ=0
v=len(p)-1-2
#for z in range(len(p)):
#    summ=summ+(((y[z]-prob[z])/error[z])**2)
chisq=(1/v)*summ
print("chisq: ",chisq)
#normmj=maxprob/ymax
#print("normmj: ",normmj)
#for h in range(len(prob)):
#    prob[h]=prob[h]/normmj
print("totalprob: ",totalprob)

print("done")
print("X VALUES: ", x)
print("Y VALUES: ", y)
a=[0.01,0.1,1,10]
#x = (line.split()[0])
#print(xs)
#print(line)
fig, ax=plt.subplots()
#ax.set_xlim(10)
#plt.xscale("log")
#ax = fig.add_subplot(111)
#ax.xaxis.set_major_locator(MultipleLocator(10))
#ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f'))
#ax = plt.gca()
#xlim = x[len(x)-1]+1
#ylim = ymax*1.5
#if logy == "1":
#    ax.set_yscale('symlog')
#    ylim = ymax*3
#if logx == "1":
 #   ax.set_xscale('symlog', linthreshx=0.1)
#ax.set_xlim([0,xlim])
#x.set_ylim([0,ylim])
#plt.xlabel("Energy (GeV)")
#plt.ylabel('Particles Per Primary')
#titlein = input("Set graph title? - \n0:no\n1:yes \n :")
#if titlein == "1":
#    title = input("Enter title: ")
#    plt.title(title)
#else:
#plt.title("Gamma rays Out - 2GeV Primaries - 2cm Target")
#ax.text(x=1.0, y=1.25, s='total particles: %d'%(int(total*float(prim))), color='#334f8d')
#ax.text(x=4, y=2e-7, s='Temperature: %dMeV'%((int(T)*k/(1.6*(10**(-19))))*10**(-6)), color='#334f8d')
#ax.text(x=1, y=ymax*1.1, s='total particles: %d'%(int(total*float(prim))), color='#334f8d')
print("this is working")
#plt.plot(x, y, label="total particles: %d"%(int(total)))
#ax.plot(x, y, 'o', color='black', markersize=2)
#(_, caps, _) = plt.errorbar(x, y, yerr=error, fmt='o',color='black', markersize=1, capsize=5)
#for cap in caps:
#    cap.set_markeredgewidth(1)
#ax.plot(x,prob, '--', color='blue', markersize=1)
#ax.xaxis.set_ticks(a)
#locmaj = matplotlib.ticker.LogLocator(base=10,numticks=5)
#ax.xaxis.set_major_locator(locmaj)
#locmin = matplotlib.ticker.LogLocator(base=10,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=9)
#ax.xaxis.set_minor_locator(locmin)
#ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
#ax.text(x=6, y=20, s="total particles: %fcm" %(total), color='#334f8d')
#plt.savefig('gamma.png', dpi=1000)
#plt.savefig('muonplot_%s_per_prim.pdf'%(search_phrase))
#plt.show()
print("total: ", total)
print(ymax)
print("Optimal T: ", T,"K")

k = 1.38064852*10**(-23)
c=3*10**8
m=1.883531627*10**(-28)
#i = input("-- specify particle type --\nenter number-\nmuons: 26\npositrons: 27\nelecrons: 28\nphotons: 29\n:" )
#search_phrase = input("-- specify data set --" )
i=28
search_phrase ="elec_Ol"
#q = open('muon_production_%d_sum.lis'%(int(i)),'r')
#line1 = q.readlines()[4]
#parts = line1.split(' ')    
prim = 1000000000
#print(prim)
#prim = input("enter number of primaries:")
total = 0
#logx = input("log plot x axis?\n0:no\n1:yes\n:")
#logy = input("log plot y axis?\n0:no\n1:yes\n:")
#prim = 93000000
error=[]

f = open('%smuon_production_%d_tab.lis'%(path,int(i)),'r')
x=[]
prob=[]
gamma=[]
K2=[]
p=[]
xs=[]
y=[]
ye=[]
xe=[]
ys=[]
errors=[]
ymax=0
norm=0.1
normmj=720962.4432
maxprob=0
print(xs)
line_num = 0
line_num_2 = 0
end_line = 100000000000000000000000000
start_line = 1000000000000
for line in f.readlines():
    index=0
    end_count = 0
    if line.find(search_phrase) >= 0:
        start_line=line_num+1
       # print(line_num)
        #print(line)
   # print("line_num:",line_num)
   # print("start_line:",start_line)
    if line_num > start_line:
         if line.find(" ", 2,3) >= 0:  
            end_line=line_num
            #print("end line:", end_line)
    if line_num > start_line and line_num < end_line:
        #print(line)
        check = (line.split())
        if len(check) == 0:
            break
        ys = float(line.split()[2])
        errors = float(line.split()[3])
        if ys != 0:
            if ymax < ys*norm:
                ymax = ys*norm
            y.append(ys*norm)
            ye.append(ys*norm)
            errore.append((ys*norm/100)*errors)
            total=total+(ys*norm)
            xs = float(line.split()[0])
            
            x.append(xs)
            xe.append(xs)
            #print("xs: ",xs)
            if ((((1.6*10**(-19))*1000000000*xs)**2)-((m*(c**2))**2))>0:
                p.append((math.sqrt((((1.6*10**(-19))*1000000000*xs)**2)-((m*(c**2))**2)))/c)
            if ((((1.6*10**(-19))*1000000000*xs)**2)-((m*(c**2))**2))<0:
                p.append(0)
            #print("p: ",p)
            K2.append(0)
            prob.append(0)
            gamma.append(0)
        #print(x)
    if line_num >= end_line:
        break
    line_num += 1
totalprob=0
k = 1.38064852*10**(-23)
c=3*10**8
m=1.883531627*10**(-28)
lowchisq=1000000000000000000000000000000
#T=[]
#for s in range(0,600,50):
#    T.append(s)
testprob=[]
testy=[]
testerror=[]
for z in range(len(p)):
    testprob.append(0)
    testy.append(0)
    testerror.append(0)
for T in range(8000000000000,8830000000000,1000000000):
    testprob=[]
    testy=[]
    testerror=[]
    for z in range(len(p)):
        testprob.append(0)
        testy.append(0)
        testerror.append(0)
    totalprob=0
    maxprob=0
    print("TEMPERATURE",T)
    theta = k*T/(m*(c**2))
    #print("theta: ",theta)
    for z in range(len(p)):
        gamma[z]=math.sqrt(1+(p[z]/(m*c))**2)
        #print("gamma: ",gamma)
        K2[z] = kn(2,1/theta)
        #K2[z]=((p[z]-1)/((gamma[0]**(-p[z]+1))-(gamma[len(gamma)-1]**(-p[z]+1))))*(gamma[z]**(-p[z]))
        #print("K2: ",K2)
        prob[z]=((((gamma[z]**2)*(math.sqrt(1-(1/(gamma[z]**2)))))/(theta*K2[z]))*math.exp(-gamma[z]/theta))
        totalprob=totalprob+prob[z]
        if prob[z]>maxprob:
            maxprob=prob[z]
        #print("prob: ",prob)
    for z in range(len(p)):
        testprob[z]=prob[z]/totalprob
        testy[z]=y[z]/total
        testerror[z]=errore[z]/total
    summ=0
    v=len(p)-1-2
    for z in range(len(p)):
        summ=summ+(((testy[z]-testprob[z])/testerror[z])**2)
    chisq=(1/v)*summ
    print("chisq: ",chisq)
    if chisq<lowchisq:
        lowchisq=chisq
        optT=T
        print("optimal T: ",optT)
T=optT
totalprob=0
maxprob=0
theta = k*T/(m*(c**2))
print("theta: ",theta)
for z in range(len(p)):
    gamma[z]=math.sqrt(1+(p[z]/(m*c))**2)
    print("gamma: ",gamma)
    K2[z] = kn(2,1/theta)
    #K2[z]=((p[z]-1)/((gamma[0]**(-p[z]+1))-(gamma[len(gamma)-1]**(-p[z]+1))))*(gamma[z]**(-p[z]))
    print("K2: ",K2)
    prob[z]=((((gamma[z]**2)*(math.sqrt(1-(1/(gamma[z]**2)))))/(theta*K2[z]))*math.exp(-gamma[z]/theta))
    totalprob=totalprob+prob[z]
    if prob[z]>maxprob:
        maxprob=prob[z]
    print("prob: ",prob)
for z in range(len(p)):
    prob[z]=(prob[z]/totalprob)*total
    y[z]=y[z]
#    error[z]=error[z]
summ=0
v=len(p)-1-2
#for z in range(len(p)):
#    summ=summ+(((y[z]-prob[z])/error[z])**2)
chisq=(1/v)*summ
print("chisq: ",chisq)
#normmj=maxprob/ymax
#print("normmj: ",normmj)
#for h in range(len(prob)):
#    prob[h]=prob[h]/normmj
print("totalprob: ",totalprob)

print("done")
print("X VALUES: ", x)
print("Y VALUES: ", y)
a=[0.01,0.1,1,10]
#x = (line.split()[0])
#print(xs)
#print(line)
fig, ax=plt.subplots()
#ax.set_xlim(10)
#plt.xscale("log")
#ax = fig.add_subplot(111)
#ax.xaxis.set_major_locator(MultipleLocator(10))
#ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f'))
#xlim = x[len(x)-1]+1
#ax = plt.gca()
#ylim = ymax*1.5
#if logy == "1":
#    ax.set_yscale('symlog')
#    ylim = ymax*3
#if logx == "1":
 #   ax.set_xscale('symlog', linthreshx=0.1)
#ax.set_xlim([0,xlim])
#x.set_ylim([0,ylim])
#plt.xlabel("Energy (GeV)")
#plt.ylabel('Particles Per Primary')
#titlein = input("Set graph title? - \n0:no\n1:yes \n :")
#if titlein == "1":
#    title = input("Enter title: ")
#    plt.title(title)
#else:
#plt.title("Gamma rays Out - 2GeV Primaries - 2cm Target")
#ax.text(x=1.0, y=1.25, s='total particles: %d'%(int(total*float(prim))), color='#334f8d')
#ax.text(x=4, y=2e-7, s='Temperature: %dMeV'%((int(T)*k/(1.6*(10**(-19))))*10**(-6)), color='#334f8d')
#ax.text(x=1, y=ymax*1.1, s='total particles: %d'%(int(total*float(prim))), color='#334f8d')
print("this is working")
#plt.plot(x, y, label="total particles: %d"%(int(total)))
#ax.plot(x, y, 'o', color='black', markersize=2)
#(_, caps, _) = plt.errorbar(x, y, yerr=error, fmt='o',color='black', markersize=1, capsize=5)
#for cap in caps:
#    cap.set_markeredgewidth(1)
#ax.plot(x,prob, '--', color='blue', markersize=1)
#ax.xaxis.set_ticks(a)
#locmaj = matplotlib.ticker.LogLocator(base=10,numticks=5)
#ax.xaxis.set_major_locator(locmaj)
#locmin = matplotlib.ticker.LogLocator(base=10,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=9)
#ax.xaxis.set_minor_locator(locmin)
#ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
#ax.text(x=6, y=20, s="total particles: %fcm" %(total), color='#334f8d')
#plt.savefig('gamma.png', dpi=1000)
#plt.savefig('muonplot_%s_per_prim.pdf'%(search_phrase))
#plt.show()
print("total: ", total)
print(ymax)
print("Optimal T: ", T,"K")

k = 1.38064852*10**(-23)
c=3*10**8
m=1.883531627*10**(-28)
#i = input("-- specify particle type --\nenter number-\nmuons: 26\npositrons: 27\nelecrons: 28\nphotons: 29\n:" )
#search_phrase = input("-- specify data set --" )
i=28
search_phrase ="posi_Ol"
#q = open('muon_production_%d_sum.lis'%(int(i)),'r')
#line1 = q.readlines()[4]
#parts = line1.split(' ')    
prim = 1000000000
#print(prim)
#prim = input("enter number of primaries:")
total = 0
#logx = input("log plot x axis?\n0:no\n1:yes\n:")
#logy = input("log plot y axis?\n0:no\n1:yes\n:")
#prim = 93000000
error=[]
f = open('%smuon_production_%d_tab.lis'%(path,int(i)),'r')
x=[]
prob=[]
gamma=[]
K2=[]
p=[]
xs=[]
xp=[]
yp=[]
y=[]
ys=[]
errors=[]
ymax=0
norm=0.1
normmj=720962.4432
maxprob=0
print(xs)
line_num = 0
line_num_2 = 0
end_line = 100000000000000000000000000
start_line = 1000000000000
for line in f.readlines():
    index=0
    end_count = 0
    if line.find(search_phrase) >= 0:
        start_line=line_num+1
       # print(line_num)
        #print(line)
   # print("line_num:",line_num)
   # print("start_line:",start_line)
    if line_num > start_line:
         if line.find(" ", 2,3) >= 0:  
            end_line=line_num
            #print("end line:", end_line)
    if line_num > start_line and line_num < end_line:
        #print(line)
        check = (line.split())
        if len(check) == 0:
            break
        ys = float(line.split()[2])
        errors = float(line.split()[3])
        if ys != 0:
            if ymax < ys*norm:
                ymax = ys*norm
            y.append(ys*norm)
            yp.append(ys*norm)
            errorp.append((ys*norm/100)*errors)
            total=total+(ys*norm)
            xs = float(line.split()[0])
            
            x.append(xs)
            xp.append(xs)
            #print("xs: ",xs)
            if ((((1.6*10**(-19))*1000000000*xs)**2)-((m*(c**2))**2))>0:
                p.append((math.sqrt((((1.6*10**(-19))*1000000000*xs)**2)-((m*(c**2))**2)))/c)
            if ((((1.6*10**(-19))*1000000000*xs)**2)-((m*(c**2))**2))<0:
                p.append(0)
            #print("p: ",p)
            K2.append(0)
            prob.append(0)
            gamma.append(0)
        #print(x)
    if line_num >= end_line:
        break
    line_num += 1
totalprob=0
k = 1.38064852*10**(-23)
c=3*10**8
m=1.883531627*10**(-28)
lowchisq=1000000000000000000000000000000
#T=[]
#for s in range(0,600,50):
#    T.append(s)
testprob=[]
testy=[]
testerror=[]
for z in range(len(p)):
    testprob.append(0)
    testy.append(0)
    testerror.append(0)
for T in range(8000000000000,8830000000000,1000000000):
    testprob=[]
    testy=[]
    testerror=[]
    for z in range(len(p)):
        testprob.append(0)
        testy.append(0)
        testerror.append(0)
    totalprob=0
    maxprob=0
    print("TEMPERATURE",T)
    theta = k*T/(m*(c**2))
    #print("theta: ",theta)
    for z in range(len(p)):
        gamma[z]=math.sqrt(1+(p[z]/(m*c))**2)
        #print("gamma: ",gamma)
        K2[z] = kn(2,1/theta)
        #K2[z]=((p[z]-1)/((gamma[0]**(-p[z]+1))-(gamma[len(gamma)-1]**(-p[z]+1))))*(gamma[z]**(-p[z]))
        #print("K2: ",K2)
        prob[z]=((((gamma[z]**2)*(math.sqrt(1-(1/(gamma[z]**2)))))/(theta*K2[z]))*math.exp(-gamma[z]/theta))
        totalprob=totalprob+prob[z]
        if prob[z]>maxprob:
            maxprob=prob[z]
        #print("prob: ",prob)
    for z in range(len(p)):
        testprob[z]=prob[z]/totalprob
        testy[z]=y[z]/total
        testerror[z]=errorp[z]/total
    summ=0
    v=len(p)-1-2
    for z in range(len(p)):
        summ=summ+(((testy[z]-testprob[z])/testerror[z])**2)
    chisq=(1/v)*summ
    print("chisq: ",chisq)
    if chisq<lowchisq:
        lowchisq=chisq
        optT=T
        print("optimal T: ",optT)
T=optT
totalprob=0
maxprob=0
theta = k*T/(m*(c**2))
print("theta: ",theta)
for z in range(len(p)):
    gamma[z]=math.sqrt(1+(p[z]/(m*c))**2)
    print("gamma: ",gamma)
    K2[z] = kn(2,1/theta)
    #K2[z]=((p[z]-1)/((gamma[0]**(-p[z]+1))-(gamma[len(gamma)-1]**(-p[z]+1))))*(gamma[z]**(-p[z]))
    print("K2: ",K2)
    prob[z]=((((gamma[z]**2)*(math.sqrt(1-(1/(gamma[z]**2)))))/(theta*K2[z]))*math.exp(-gamma[z]/theta))
    totalprob=totalprob+prob[z]
    if prob[z]>maxprob:
        maxprob=prob[z]
    print("prob: ",prob)
for z in range(len(p)):
    prob[z]=(prob[z]/totalprob)*total
    y[z]=y[z]
#    error[z]=error[z]
summ=0
v=len(p)-1-2
#for z in range(len(p)):
#    summ=summ+(((y[z]-prob[z])/error[z])**2)
chisq=(1/v)*summ
print("chisq: ",chisq)
#normmj=maxprob/ymax
#print("normmj: ",normmj)
#for h in range(len(prob)):
#    prob[h]=prob[h]/normmj
print("totalprob: ",totalprob)

print("done")
print("X VALUES: ", x)
print("Y VALUES: ", y)
a=[0.01,0.1,1,10]
#x = (line.split()[0])
#print(xs)
#print(line)
fig, ax=plt.subplots()
#ax.set_xlim(10)
#plt.xscale("log")
#ax = fig.add_subplot(111)
#ax.xaxis.set_major_locator(MultipleLocator(10))
#ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f'))
ax = plt.gca()
#xlim = x[len(x)-1]+1
#ylim = ymax*1.5
#if logy == "1":
#    ax.set_yscale('symlog')
#    ylim = ymax*3
#if logx == "1":
 #   ax.set_xscale('symlog', linthreshx=0.1)
#ax.set_xlim([0,xlim])
#x.set_ylim([0,ylim])

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)

ax.set_yscale('log')
plt.xlabel("Energy (GeV)")
plt.ylabel('Particles Per Primary')
#titlein = input("Set graph title? - \n0:no\n1:yes \n :")
#if titlein == "1":
#    title = input("Enter title: ")
#    plt.title(title)
#else:
#plt.title("Gamma rays Out - 2GeV Primaries - 2cm Target")
#ax.text(x=1.0, y=1.25, s='total particles: %d'%(int(total*float(prim))), color='#334f8d')
#ax.text(x=4, y=2e-7, s='Temperature: %dMeV'%((int(T)*k/(1.6*(10**(-19))))*10**(-6)), color='#334f8d')
#ax.text(x=1, y=ymax*1.1, s='total particles: %d'%(int(total*float(prim))), color='#334f8d')
print("this is working")
ax.legend()
#plt.plot(x, y, label="total particles: %d"%(int(total)))
ax.plot(xp, yp, 'o', color='red', markersize=2, marker="o", label="positron")
#(_, caps, _) = plt.errorbar(xp, yp, yerr=errorp, fmt='o',color='black', markersize=1, capsize=5)
#for cap in caps:
#    cap.set_markeredgewidth(1)
ax.plot(xe, ye, 'o', color='blue', markersize=2, marker = "^", label="electron")
#(_, caps, _) = plt.errorbar(xe, ye, yerr=errore, fmt='o',color='black', markersize=1, capsize=5)
#for cap in caps:
#    cap.set_markeredgewidth(1)
ax.plot(xg, yg, 'o', color='green', markersize=2, marker="+", label="gamma")
#(_, caps, _) = plt.errorbar(xg, yg, yerr=errorg, fmt='o',color='black', markersize=1, capsize=5)
#for cap in caps:
#    cap.set_markeredgewidth(1)
leg = ax.legend();
matplotlib.rc('xtick') 
matplotlib.rc('ytick')
plt.rc('legend',fontsize=13) 

#(_, caps, _) = plt.errorbar(x, y, yerr=error, fmt='o',color='black', markersize=1, capsize=5)
#for cap in caps:
#    cap.set_markeredgewidth(1)
#ax.plot(x,prob, '--', color='blue', markersize=1)
#ax.xaxis.set_ticks(a)
#locmaj = matplotlib.ticker.LogLocator(base=10,numticks=5)
#ax.xaxis.set_major_locator(locmaj)
#locmin = matplotlib.ticker.LogLocator(base=10,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=9)
#ax.xaxis.set_minor_locator(locmin)
#ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
#ax.text(x=6, y=20, s="total particles: %fcm" %(total), color='#334f8d')
plt.savefig('%snoise_out.png'%(path), bbox_inches='tight', dpi=1000)
#plt.savefig('muonplot_%s_per_prim.pdf'%(search_phrase))
plt.show()
#print("total: ", total)
#print(ymax)
#print("Optimal T: ", T,"K")
print(xg)
print(yg)
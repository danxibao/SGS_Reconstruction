
import os
import random as rd
import math


if not os.path.isdir("Detection"):
    os.mkdir("Detection")
    
f1=open("Detection\\_SourceDefine.txt",'w')
f2=open("Detection\\0_RamdomData.txt",'w')
f3=open("Detection\\0_Activity.txt",'w')



radius=0#cm
angle=0#degree
height=0#cm
NoEnergy=[0,1,2,5]#MeV
activity=10000000#bq
SumActivity=0

activitytemp=0

#设置
#废物桶半径，cm
ROfWD=35
#废物桶高度，cm
HOfWD=105

#发射源是否随机分布
IsRandom=True
#密度，g/cm3
density=0.5
#测量方式//SGS ISGS STGS2ER STGS2EB STGS4ER STGS4EB STGS8ER STGS8EB *****
DetectionType="SGS"
#单点源分布时点源高度
height=47.25

#多点源随机分布99个废物桶，单点源分布14个废物桶
NuDrums=99 if IsRandom else 14
#多点源随机分布时设置每个桶内点源数
NuSources=35
#输出
for nDrum in range(NuDrums):
    f1.write("$$Drum"+str(nDrum+1)+"\n")
    f1.write("Detection:"+str(DetectionType)+";\n")
    f1.write("Density:"+str(density)+";\n")
    SumActivity=0
    
    if not IsRandom:
        NuSources=4
        radius=1.0*nDrum*ROfWD/NuDrums
        for nSource in range(NuSources):
            f1.write("##Source"+str(nSource+1)+'\n')
            f1.write("Radius:"+str(radius)+";\n")
            f1.write("Angle:"+str(angle)+";\n")
            f1.write("Height:"+str(height)+";\n")
            f1.write("NoEnergy:"+str(NoEnergy[nSource])+";\n")
            f1.write("Activity:"+str(activity)+";\n")
            SumActivity+=activity
    else:
        for nSource in range(NuSources):
            if nSource%4==0:
                
                tmp=rd.randint(0,999)
                #rd.seed(10)
                rand1=HOfWD*((rd.randint(0,9999)*tmp)%10000)/10000.0
                #rd.seed(121)
                rand2=ROfWD*round(math.sqrt(((rd.randint(0,9999)*tmp)%10000)/10000.0),4)
                #rd.seed(521)
                rand3=((rd.randint(0,9999)*tmp)%10000)/10000.0
                height=rand1
                radius=rand2
                activitytemp=activity*(0.25+0.75*rand3)
                SumActivity+=activitytemp
                f2.write(str(rand1)+" "+str(rand2)+" "+str(rand3)+"\n")
            
            f1.write("##Source"+str(nSource+1)+'\n')
            f1.write("Radius:"+str(radius)+";\n")
            f1.write("Angle:"+str(angle)+";\n")
            f1.write("Height:"+str(height)+";\n")
            f1.write("NoEnergy:"+str(NoEnergy[nSource%4])+";\n")
            f1.write("Activity:"+str(activitytemp)+";\n")
    
    f1.write("**SumA:"+str(SumActivity)+";\n")
    f1.write("\n")
    f3.write(str(SumActivity)+"\n")
f1.write("End\n")
    
                    
                       
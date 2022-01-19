# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 00:43:06 2019

@author: batuhan
"""

#section with multi layer reinforcements.
# Input variables
import numpy as np
import pandas as pd
import read_excelData as readE
from pandas import DataFrame

#from read_excelData import read_excel


def calc_momentCapacity():
    
    n_ex,list_beam = readE.read_excel()
  
    
    
    BAR_SIZE_CONSTANT = 2
    
    # GEOMETRY PROPERTIES :
    h=60             # h: height of the RC beam(in cm)
    bw=30            # bw: width of the RC beam(in cm)               
    
    clear_cover=4    # clear_cover is the distance btw lower face of the (in cm)
                      # bar to the concrete surface.
    d_eff=h-clear_cover   
    Mr_span=np.zeros((n_ex))
    A=np.zeros((4,1))  
    
    As_req=np.zeros(n_ex,dtype=float)
    Md=np.zeros(n_ex,dtype=float)
    
    
    Md_left=np.zeros(n_ex,dtype=float)
    As_left_req=np.zeros(n_ex,dtype=float)  
    As_left=np.zeros(n_ex,dtype=float)  
    Mr_left=np.zeros(n_ex,dtype=float)
    
    As_left_tension=[]
    As_left_compression=[]
    As_right_tension=[]
    As_right_compression=[]
    
    Fİ=[]
    As_REQ=[]
    Md_SPAN=[]
    Mr_SPAN=[]
    beam_name=[]
    Md_right=np.zeros(n_ex,dtype=float)
    As_right_req=np.zeros(n_ex,dtype=float)  
    As_right=np.zeros(n_ex,dtype=float)  
    Mr_right=np.zeros(n_ex,dtype=float)
    
    Mr_supports=np.zeros((n_ex,4),dtype=float)
    c_c=np.zeros(n_ex,dtype=float)
    C=[]
    bar_size=12
    bool_controller = False
    for i in range(0,n_ex,1):
        fi_n_d=np.array([[bar_size,2,clear_cover],   # first column: diameter of steel bar(in cm) exp.fi20
                 [bar_size,2,h/2],              # second column: number of steel bar
                 [0,0,0],              # third column: d(i)(effective depths)            
                 [bar_size,3,h-clear_cover]])
        
        fi=((fi_n_d[:,0])*(10**-1))   # diameter of the steel bars(in cm) 
        n=fi_n_d[:,1]      # n is the steel bars number for each layer.    
        nx=len(n)
        
        for m in range(0,4,1):
            
            As=(((np.pi)*(n[m])*((fi[m])**2))/4)   # As: Cross sectional area of the steel bars
            A[m] = As
            sizeA=len(A)
      
    
        ro=A[3]/(h*bw) # Reinforcement Tension bars ratio
        
    
        # MATERIAL PROPERTIES :
    
        conc_class= "C25"
        string_fck=conc_class[1:] 
        fck=int(string_fck)       # Characteristic Strength of Concrete(in ton/cm^2)      
        fcd=((fck*(10**-2))/1.5)  # Design Strength of Concrete(in ton/cm^2)   
        steel_class="S420"
        string_fyk=steel_class[1:]
        fyk=int(string_fyk)       # Characteristic Strength of Steel (in ton/cm^2)       
        fyd=((fyk*(10**-2))/1.15) # Design Strength of Steel (in ton/cm^2)
    
        Es=2*(10**6)              # Modulus of Elasticity of Steel (in kgf/cm2)  
        fctk=0.35*np.sqrt(fck)
        fctd=fctk/1.5
        ro_min=0.8*(fctd/(fyk/1.15))
    
        # Check the steel ratio with respect to the TS500
        if ro>ro_min and ro<=0.02:
            print(">>> (ro) steel ratio higher than min steel ratio(ro_min)","ro=",ro,">","ro_min=",ro_min)         
            print(">>> (ro) steel ratio smaller than 0.02","ro=",ro,"<","0.02")
              
        st_balanced_ratio=([["S220","C16","0.0316"],     # The values come from Ersyo's book.
                        ["S220","C18","0.0344"],     # First column is steel class
                        ["S220","C20","0.0373"],     # Second column is concrete class
                        ["S220","C25","0.0488"],     # Third column is balanced ratios
                        ["S420","C16","0.0135"],
                        ["S420","C20","0.0148"],
                        ["S420","C25","0.0209"],
                        ["S420","C26","0.0237"],
                        ["S420","C35","0.0263"],
                        ["S420","C40","0.0297"],
                        ["S420","C45","0.0317"],
                        ["S420","C50","0.0334"]])
    
        # Balanced ratio is taken from s_balanced _ratio.
        size1=len(st_balanced_ratio)
        balanced_ratio_string="invalid!"
        for b in range(0,size1,1):
            if (steel_class==st_balanced_ratio[b][0]) and (conc_class==st_balanced_ratio[b][1]):
                balanced_ratio_string=st_balanced_ratio[b][2]
        if balanced_ratio_string=="invalid!":
            print("Invalid Steel and Concrete Class!")
        else:
            ro_b=float(balanced_ratio_string)
        
        
        
        Md[i]=list_beam[i][2]
        Md_SPAN.append(Md[i])
        As_req[i]=(Md[i]*10)/(fyd*0.7*(h-clear_cover))   #cm2
        if As_req[i]<0:
            As_req[i]=As_req[i]*(-1)
        else:
            As_req[i]=As_req[i]*(1)
        
        if Md[i]<0:
            Md[i]=-1*Md[i]
        else:
            Md[i]=Md[i]
        
        As_REQ.append(As_req[i])
        if A[3]>As_req[i]:
            if(bool_controller):
                bar_size-=BAR_SIZE_CONSTANT
                bool_controller=False
 
                
            Fİ.append(bar_size)
                
            #print("Moment Capacity is OK for beam number", int(list_beam[i][0]),"Mr=",M*10**-1,">","Md=",Md[i])
            print("Area of tension bars are higher than required area","A=",A[3],"cm2",">","As_req=",As_req[i])
            
        else:
            print("INCREASE THE BAR SIZE OR INCREASE THE NUMBER OF BAR")
            
            bar_size=bar_size+BAR_SIZE_CONSTANT
            
            bool_controller=True
            
            Fİ.append(bar_size)
          
        fi_n_d=np.array([[bar_size,2,clear_cover],   # first column: diameter of steel bar(in cm) exp.fi20
                 [bar_size,2,h/2],              # second column: number of steel bar
                 [0,0,0],              # third column: d(i)(effective depths)            
                 [bar_size,3,h-clear_cover]])
        
        fi=((fi_n_d[:,0])*(10**-1))   # diameter of the steel bars(in cm) 
        n=fi_n_d[:,1]      # n is the steel bars number for each layer.    
        nx=len(n)
        
        for m in range(0,4,1):
            
            As=(((np.pi)*(n[m])*((fi[m])**2))/4)   # As: Cross sectional area of the steel bars
            A[m] = As
            sizeA=len(A)
        
        As_left_tension.append(A[3][0]) 
        As_left_compression.append(A[0][0])
        
        As_right_tension.append(A[3][0])
        As_right_compression.append(A[0][0])
        #For Stirrup Calculation 
        
        # k1 value is taken from conc_class_k1    
        conc_class_k1=([["C16","0.85"],   # First column is concrete class
                    ["C18","0.85"],   #Second column is k1 values.
                    ["C20","0.85"],
                    ["C25","0.85"],
                    ["C26","0.82"],
                    ["C35","0.79"],
                    ["C40","0.76"],
                    ["C45","0.73"],
                    ["C50","0.70"]])
        size=len(conc_class_k1)
        k1_string="invalid"
        for s in range(0,size,1):
            if conc_class==conc_class_k1[s][0]:
                k1_string=conc_class_k1[s][1]
        if k1_string=="invalid":
            print("Invalid Concrete Class!")
        else:
            k1=float(k1_string)
        
    
            d=fi_n_d[:,2]        
            dx=len(d)
            epsilon_sy=((fyd*(10**3))/Es)
            fs_matrix=np.zeros((dx,1))
        
        
            start_c=clear_cover+1    # Assuming first c value is clear_cover+1        
            controller=0             # controller controls the range of c to arrange the equation.
        
            for j in np.arange(start_c,d[dx-1],0.01): 
                c=j
                
                if d[0]<c and c<d[1]:
                    epsilon1=(0.003*((c-d[0])/(c)))
               
        
                    if (epsilon1>epsilon_sy):
                        fs1=fyd
                    else:
                        fs1=(epsilon1*Es)*(10**-3)
                    fs_matrix[0]=fs1
                    for k in range(1,dx,1):
                        epsilon=(0.003*((d[k]-c)/(c)))
                        if (epsilon>epsilon_sy):
                            fs=fyd
                        else:
                            fs=(epsilon*Es)*(10**-3)
                        fs_matrix[k]=fs
                        controller=1
                  
                elif (d[1]<c and c<d[2]):
                
                    for p in range(0,2,1):
                        epsilon=(0.003*((c-d[p])/(c)))
                        if (epsilon>epsilon_sy):
                            fs=fyd
                        else:
                            fs=(epsilon*Es)*(10**-3)
                        fs_matrix[p]=fs
                    for r in range(dx-2,dx,1):
                        epsilon=(0.003*((d[r]-c)/(c)))
                        if (epsilon>epsilon_sy):
                            fs=fyd
                        else:
                            fs=(epsilon*Es)*(10**-3)
                        fs_matrix[r]=fs
                        controller=2
            
                elif (d[2]<c and c<d[3]):
                    for a in range(0,dx-1,1):
                        epsilon=(0.003*((c-d[a])/(c)))
                        if (epsilon>epsilon_sy):
                            fs=fyd
                        else:
                            fs=(epsilon*Es)*(10**-3)
                        fs_matrix[a]=fs
                        
                    d4=d[a+1]
                    epsilon=(0.003*((d4-c)/(c)))
                    if (epsilon>epsilon_sy):
                        fs=fyd
                    else:
                        fs=(epsilon*Es)*(10**-3)
                    fs_matrix[dx-1]=fs
                    controller=3
    
            
                Fs1=fs_matrix[0]*A[0]
                Fs2=fs_matrix[1]*A[1]
                Fs3=fs_matrix[2]*A[2]
                Fs4=fs_matrix[3]*A[3]
            
                if controller==1:
                    sum_Fx=(0.85*fcd*k1*c*bw)+Fs1-Fs2-Fs3-Fs4
                elif controller==2:
                    sum_Fx=(0.85*fcd*k1*c*bw)+Fs1+Fs2-Fs3-Fs4
                elif controller==3:
                    sum_Fx=(0.85*fcd*k1*c*bw)+Fs1+Fs2+Fs3-Fs4
                else:
                    break
            
                if 0<sum_Fx and sum_Fx<1:
                    break
            #print(">>> c value is equal to:",round(c,2),"cm")    
            c_c[i]=round(c,2)
          
            if controller==1:
                if ((k1*c)/2)<d[0]:
                    Mr_span[i]=float(Fs1*(d[0]-((k1*c)/2)))+(Fs2*(d[1]-((k1*c)/2)))+(Fs3*(d[2]-((k1*c)/2)))+(Fs4*(d[3]-((k1*c)/2))) # casting
                    Mr_supports[i][3]=(Mr_span[i])*10**-1
                else:
                    Mr_span[i]=float((Fs1*(((k1*c)/2)-d[0]))+(Fs2*(d[1]-((k1*c)/2)))+(Fs3*(d[2]-((k1*c)/2)))+(Fs4*(d[3]-((k1*c)/2))))
                    Mr_supports[i][3]=(Mr_span[i])*10**-1
            elif controller==2:
                if ((k1*c)/2)<d[1]:
                    Mr_span[i]=float(Fs1*(((k1*c)/2)-d[0]))+(Fs2*(d[1]-((k1*c)/2)))+(Fs3*(d[2]-((k1*c)/2)))+(Fs4*(d[3]-((k1*c)/2)))
                    Mr_supports[i][3]=(Mr_span[i])*10**-1
                else:
                    Mr_span[i]=float(Fs1*(((k1*c)/2)-d[0]))+(Fs2*(((k1*c)/2)-d[1]))+(Fs3*(d[2]-((k1*c)/2)))+(Fs4*(d[3]-((k1*c)/2)))
                    Mr_supports[i][3]=(Mr_span[i])*10**-1
            elif controller==3:
                if ((k1*c)/2)<d[2]:
                    Mr_span[i]=float(Fs1*(((k1*c)/2)-d[0]))+(Fs2*(((k1*c)/2)-d[1]))+(Fs3*(d[2]-((k1*c)/2)))+(Fs4*(d[3]-((k1*c)/2)))
                    Mr_supports[i][3]=(Mr_span[i])*10**-1
                else:
                    Mr_span[i]=float(Fs1*(((k1*c)/2)-d[0]))+(Fs2*(((k1*c)/2)-d[1]))+(Fs3*(((k1*c)/2)-d[2]))+(Fs4*(d[3]-((k1*c)/2)))
                    Mr_supports[i][3]=(Mr_span[i])*10**-1
            #print(">>> Moment Capacity is: ",Mr_span[i]*10*10**-2,"kN.m")
            
            #Balanced reinforcement check
                
            if fs_matrix[0]==fyd:
                check_value=((A[sizeA-1]-A[0])/(bw*d[dx-1]))
            else:
                check_value=(((A[sizeA-1])/(bw*d[dx-1]))-(((A[0])/(bw*d[dx-1]))*((fs_matrix[0])/fyd)))
        
        
            if (check_value<ro_b):
                print(">>> Balanced Reinforcement Check is OK")
            else:
                print(">>> Balanced Reinforcement Check is not OK")
    
            if Mr_span[i]>Md[i]:
                print("Moment Capacity is OK for beam number", int(list_beam[i][0]),"Mr_span=",Mr_span[i]*10**-1,">","Md=",Md[i])
            #print("Area of tension bars are higher than required area","A=",A[3],"cm2",">","As_req=",As_req)
                
            else:
                print("INCREASE THE BAR SIZE OR INCREASE THE NUMBER OF BAR")
             
# Calculation of additional support reinforcement bars
        

        fi_sup=np.array([[1.2],[1.4],[1.6],[1.8],[2.0],[2.2],[2.4],[2.6],[2.8]]) 
    
        n_sup=len(fi_sup)
        
    
        Md_left[i]=list_beam[i][1]
        Md_right[i]=list_beam[i][3]
        
        if Md_left[i]<0:
            Md_left[i]=(-1)*Md_left[i]
  
        else:
            Md_left[i]=1*Md_left[i]
            
            
        if Md_right[i]<0:
            Md_right[i]=(-1)*Md_right[i]
        else:
            Md_right[i]=1*Md_right[i]
        
        
        a_left=d_eff-(np.sqrt((d_eff**2)-(2*(Md_left[i]*(10**2)))/(0.85*(fcd*10)*bw)))    
        As_left_req[i]=Md_left[i]*10**2/(fyd*10*(d_eff-(a_left/2)))

        a_right=d_eff-(np.sqrt((d_eff**2)-(2*(Md_right[i]*(10**2)))/(0.85*(fcd*10)*bw)))    
        As_right_req[i]=Md_right[i]*10**2/(fyd*10*(d_eff-(a_right/2)))


        for x in range(0,n_sup,1):
            As_left[i]=(((np.pi)*((fi_sup[x][0])**2))/4)*2
            
            
            if As_left[i]>As_left_req[i]:
                print("As_left > As_left_req",As_left[i],">",As_left_req[i], "Using 2 number",fi_sup[x]*10,"mm bar")
                break
            
            As_right[i]=(((np.pi)*((fi_sup[x][0])**2))/4)*2
            if(As_right[i]>As_right_req[i]):
                print("As_right > As_right_req",As_right[i],">",As_right_req[i], "Using 2 number",fi_sup[x]*10,"mm bar")
                break
            
   
        
        totalLeft_are=A[3]+As_left[i]
        k1_c_left=((totalLeft_are)*(fyd*10))/(0.85*fcd*10*bw)
        Mr_left[i]=totalLeft_are*fyd*10*(d_eff-(k1_c_left/2))
            
        if Mr_left[i]*10**-2>Md_left[i]:
            print("Mr_left > Md_left","Mr_left=",Mr_left[i]*10**-2,"kN.m",">","Md_left=",Md_left[i],"kN.m")
                
            Mr_supports[i][1]=Mr_left[i]*10**-2
          
        totalRight_are=A[3]+As_right[i]
        k1_c_right=((totalRight_are)*(fyd*10))/(0.85*fcd*10*bw)
        Mr_right[i]=totalRight_are*fyd*10*(d_eff-(k1_c_right/2))
            
        
        if Mr_right[i]*10**-2>Md_right[i]:
            print("Mr_right > Md_right","Mr_right=",Mr_right[i]*10**-2,"kN.m",">","Md_right=",Md_right[i],"kN.m")
        
            Mr_supports[i][2]=Mr_right[i]*10**-2
            Mr_supports[i][0]= int(list_beam[i][0]) 
            beam_name.append(Mr_supports[i][0])
            Mr_SPAN.append(Mr_supports[i][3])
            #fi.append(bar_size)
            
            
            a = DataFrame({'Beam Number': beam_name,'h': h, 'bw': bw , 'clear_cover': clear_cover,'d(effective depth':h-clear_cover,
                           'Mr_SPAN': Mr_SPAN, 'Md_SPAN': Md_SPAN,'ϕ': Fİ,})
            writer = pd.ExcelWriter('BeamDesignMoments.xlsx', engine='xlsxwriter')
            #writer.set_font_name("Times New Roman")
            a.to_excel(writer, sheet_name='Sheet1')
            worksheet = writer.sheets['Sheet1']
            worksheet.write_string(0, 0, '')
        writer.save()
        

        fi_plt=np.array([[bar_size,clear_cover,clear_cover],
                         [bar_size, bw/2, clear_cover],
                         [bar_size,bw-clear_cover,clear_cover],
                 [bar_size,clear_cover,h/2],
                 [bar_size,bw-clear_cover,h/2],
                 [bar_size,clear_cover,h-clear_cover],
                 [bar_size,bw-clear_cover,h-clear_cover]])
    
        
        
        
        pp=len(fi_plt)
        
        
        print("c values for each beam:",c_c,"cm")
        
        print(Mr_supports)
  
        import matplotlib.pyplot as plt    
        from matplotlib.patches import Rectangle
    
        #define Matplotlib figure and axis
        fig, ax = plt.subplots()
        ax.set(xlim=(0, 80), ylim = (0, 80))
    
        #create simple line plot
        ax.plot([0, 0],[0, 0])
        #add rectangle to plot
        ax.add_patch(Rectangle((0, 0), bw, h,edgecolor = 'pink',
             facecolor = 'lightsteelblue',
             fill=True,
             lw=clear_cover))  
        for y in range(0,pp,1):
    
            circle = plt.Circle((fi_plt[y][1], fi_plt[y][2]), ((fi_plt[y][0])/2)*10**-1,color='black',label="bar_size" )
            ax.add_artist(circle)
            plt.title(Mr_supports[i][0])
            
        plt.show()
                
        #fig.savefig('plotcircles.png')
 
 

calc_momentCapacity()
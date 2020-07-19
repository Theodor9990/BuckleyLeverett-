"""                               Buckley Leverett equation for fractional water calculations.

Other calculations includes:
    
    - Oil produced in PV vs. water Injected
    - Average water saturation after BT vs. water injected in PV
    - Water cut vs. oil produced in PV
    - Fractional flow as a function of oil viscosity
    - Relative permeability curves segmented for OW, WW and MixedW systems.
    - Recovery factor vs. water injected. 

To generate the relative perm. curves it is used Corey model                                    """

"-------------------------------------------------------------------------------------------------"
#Getting the data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


"-------------------------------------------------------------------------------------------------"

"""                                    A LIST OF VARIABLES                                       """

"""Here only change the name in the paranthesis for example pd.read_excel("NAME.xlsx")"""
dataset = pd.read_excel("Buckley_Leverett_Python_Input.xlsx")

"""Get viscosities array"""
viscosities_name1 = dataset.iloc[11, 14:]
viscosities_name = []
for i in viscosities_name1:
    if i>=0:
        viscosities_name.append(i)
        

""" OIL - WATER SYSTEM """                     
swc = dataset.iloc[8,3]                                  #residual water saturation
timestep = dataset.iloc[12, 3]                            #timestep(preferably do not change)                               
sor_w = dataset.iloc[8,4]                                #residual oil saturation
Sw = np.linspace(swc, 1-sor_w, timestep - 1)             #water saturation


#Relative permeability curves using Corey model                           
krw_end = dataset.iloc[9,3]                              #end water point 
krow_end = dataset.iloc[9,4]                             #end oil-water point 
n_w = dataset.iloc[10, 3]                                 #Corey exponent water
n_ow = dataset.iloc[10, 4]                                #Corey exponent oil /w water
                                   
"-------------------------------------------------------------------------------------------------"
#Generating krw
krw = []
for i in Sw:
    if i < int(1):
        krw.append(krw_end * ((i - swc)/(1 - swc - sor_w))**n_w)
    else:
        krw.append(1)
krw = np.array(krw)

## Generating kro
krow = []
for i in Sw:
    if (1 - i) < sor_w:
        krow.append(0)
    else:
        krow.append(krow_end*((1 - i - sor_w)/(1 - swc - sor_w))**n_ow)
krow = np.array(krow)

"-------------------------------------------------------------------------------------------------"

""" OIL - GAS  SYSTEM """
# Data input for the oil and gas system
sgc = dataset.iloc[8,8]                                          #residual gas saturation
sor_g = dataset.iloc[8,9]                                        #residual oil saturation
Sg = np.linspace(sgc, 1 - swc - sor_g, timestep)     


#Relative permeability curves using Corey model                  
krg_end = dataset.iloc[9,8]              #end gas point 
krog_end = dataset.iloc[9,9]             #end oil-gas point 
ng = dataset.iloc[10,8]                   #Corey exponent gas
nog = dataset.iloc[10,9]                  #Corey exponent oil /w gas

#Generating relative gas peremeability curve
krg = []
for i in Sg:
    if float(i) == 0:
        krg.append(0)
    elif float(i) >= (1 - swc):
        krg.append(krg_end)
    else: 
        krg.append(krg_end*((i - sgc)/(1 - swc - sor_g - sgc))**ng)
krg = np.array(krg)

#Generating krog
krog = []
for i in Sg:
    if float(i) == 0:
        krog.append(krow_end)
    elif 1 - float(i) - swc <= sor_g:
        krog.append(0)
    else: 
        krog.append(krow_end*((1 - float(i) - swc - sor_g)/(1 - swc - sor_g))**nog)
krog = np.array(krog)

"-------------------------------------------------------------------------------------------------"
"-------------------------------------------------------------------------------------------------"

"""Generating the water front equations""" 

# Generating fw 1
# Viscosities array
miu_w1 = dataset.iloc[12, 14:].values
miu_o1 = dataset.iloc[11, 14:].values
miu_w = []
miu_o = []

for i in miu_w1:
    if i>=0:
        miu_w.append(i)
        
for i in miu_o1:
    if i>=0:
        miu_o.append(i)

fw_1 = []
for j in range(len(miu_w)):
    fw_1.append(0)
    for i in range(len(krw)):
        if krw[i] == 0 and fw_1[i-1] != 1 :
            fw_1.append(0)
        elif krow[i] == 0 or fw_1[i-1] == 1:
            fw_1.append(1)
        else:
            fw_1.append(1/(1+(miu_w[j]*krow[i])/(miu_o[j]*krw[i])))

fw_1 = np.array(fw_1)
fw_1 = np.reshape(fw_1, newshape = (len(miu_o), timestep))

"-------------------------------------------------------------------------------------------------"

""" In order to generate the fw and perform all the logic gates, it was require to
    increase the number of steps for the matrix as it always had to check with the previous value.
    Hence, a zero was added to the beginning of the array. This will not generate issues in the actual
    code, but it is required to add a zero to the Sw at the beginning in order to uniformalize the length 
    of the arrays. Hence a new Sw is defined, containing an extra integer of a value zero 
    at the beginning for plotting the Sw vs fw """
    
Sw_fw = np.insert(Sw, 0, 0)

# Calculating the derivatives of fw/Sw
miu_w_plot = dataset.iloc[21, 18] # water viscosity 
miu_o_plot = dataset.iloc[20, 18] # oil viscosity

Sw_fw_plot = np.insert(Sw_fw, 5*[0], 0)

fw_plot = [0]

for i in range(len(krw)):
    if krw[i] == 0 and fw_plot[i-1] != 1:
        fw_plot.append(0)
    elif krow[i] == 0 or fw_plot[i-1] == 1:
        fw_plot.append(1)
    else:
        fw_plot.append(1/(1+(miu_w_plot*krow[i])/(miu_o_plot*krw[i])))

fw_plot = np.array(fw_plot)
fw_plot = np.insert(fw_plot, 5*[0], 0)

"-------------------------------------------------------------------------------------------------"

# Finding the tangent. Initialising dfw / dSw 
dfw_dSw = []
for i in range(len(Sw_fw_plot)):
    if Sw_fw_plot[i] <= swc:
        dfw_dSw.append(0)
    else:
        dfw_dSw.append(fw_plot[i]/(Sw_fw_plot[i] - swc))    
dfw_dSw = np.array(dfw_dSw)
        
"-------------------------------------------------------------------------------------------------"

# Derivative of the fw function
fw_derivative = [0]
for i in range(5, len(fw_plot) - 1):
    if Sw_fw_plot[i] and Sw_fw_plot[i-1] == 0:
        fw_derivative.append(0)
    elif Sw_fw_plot[i+1] and Sw_fw_plot[i] == 0:
        fw_derivative.append(0)
    else:
        fw_derivative.append(0.5*(fw_plot[i] - fw_plot[i-1])/(Sw_fw_plot[i]-Sw_fw_plot[i-1]) 
        + 0.5*((fw_plot[i+1] - fw_plot[i])/(Sw_fw_plot[i+1]-Sw_fw_plot[i])))
    
fw_derivative = np.array(fw_derivative)
fw_derivative = np.insert(fw_derivative, 4*[0], 0)   
fw_derivative = np.insert(fw_derivative, len(fw_derivative) , 0)   

"-------------------------------------------------------------------------------------------------"

# Gradient of the derivative of the fw function
negative_dfw_dSw = []

for i in range(len(Sw_fw_plot)):
    if Sw_fw_plot[i] < swc:
        negative_dfw_dSw.append(0)
    else:
        negative_dfw_dSw.append(fw_derivative[i] - dfw_dSw[i] )

negative_dfw_dSw = np.array(negative_dfw_dSw)

"-------------------------------------------------------------------------------------------------"

# Getting the Sw at the maximum global of the fw function
Sw_at_dfw_dSw =[]
for i in range(len(negative_dfw_dSw)):
    if negative_dfw_dSw[i] > 0 and negative_dfw_dSw[i+1] < 0:
        Sw_at_dfw_dSw.append(Sw_fw_plot[i] + ((Sw_fw_plot[i+1] - Sw_fw_plot[i])/(negative_dfw_dSw[i+1] - negative_dfw_dSw[i]))*(0-negative_dfw_dSw[i]))
        
# Getting the fw at the maximum global of the fw function
fw_at_dfw_dSw =[]
for i in range(len(negative_dfw_dSw)):
    if negative_dfw_dSw[i] > 0 and negative_dfw_dSw[i+1] < 0:
        fw_at_dfw_dSw.append(fw_plot[i] + ((fw_plot[i+1] - fw_plot[i])/(negative_dfw_dSw[i+1] - negative_dfw_dSw[i]))*(0-negative_dfw_dSw[i]))
        
"-------------------------------------------------------------------------------------------------"
        
# Sw at fw_1 
# Drowing the tangent
Sw_at_fw_1 = swc + (Sw_at_dfw_dSw[0] - swc)/fw_at_dfw_dSw[0]

# Tangent array
x_tangent = [swc, Sw_at_dfw_dSw[0], Sw_at_fw_1]
y_tangent = [0, fw_at_dfw_dSw[0], 1 ]

"-------------------------------------------------------------------------------------------------"

# Average water saturation behind shock front at wbt
Sw_avg_behind_front = swc + ((Sw_at_dfw_dSw[0] - swc)/fw_at_dfw_dSw[0])

"-------------------------------------------------------------------------------------------------"

# Creating the arrays after WBT
Swe = [Sw_at_dfw_dSw[0]]
for i in Sw_fw_plot:
    if i >= Sw_at_dfw_dSw[0]:
        Swe.append(i)
  
fwe = [fw_at_dfw_dSw[0]]
for i in fw_plot:
    if i >= fw_at_dfw_dSw[0]:
        fwe.append(i)

d_Swe = []
for i in range(len(Swe) - 1):
    d_Swe.append(Swe[i+1] - Swe[i])

d_fwe = []
for i in range(len(fwe) - 1):
    d_fwe.append(fwe[i+1] - fwe[i])
      
d_fwe_d_Swe = []
for i in range(len(d_fwe)):
    d_fwe_d_Swe.append(d_fwe[i]/d_Swe[i]) 
    
Swe_av = []
for i in range(len(Swe) - 1):
    Swe_av.append((Swe[i+1] + Swe[i])/2)

wc = []
for i in range(len(fwe) - 1):
    wc.append((fwe[i+1] + fwe[i])/2)
    
"-------------------------------------------------------------------------------------------------"
    
# Water injected in fraction of pore volume PV / Wid
Wid = []
for i in d_fwe_d_Swe:
    if i != 0:
        Wid.append(1/i)
    
# Average water saturation behind shock front after Wbt
Sw_average = []

for i in range(len(Swe_av)):
    if Swe_av[i] + (1 - fwe[i])/d_fwe_d_Swe[i] <= 1 - int(sor_w):
        Sw_average.append(Swe_av[i] + (1 - fwe[i])/d_fwe_d_Swe[i])

# Oil produced in fractions of pore volume
Npd = []
for i in Sw_average:
    Npd.append(i-swc)

# Recovery factor
  
rf = []
for i in Sw_average:
    rf.append((i-swc)/(1-swc))
    
"-------------------------------------------------------------------------------------------------"

Wid_oil_wet_reference = [0, 0.1092, 0.24, 0.81, 1.4, 2.38, 2.5, 3.5]
rf_oil_wet_reference = [0, 0.1738,  0.27, 0.30, 0.33, 0.37, 0.375, 0.38]

Wid_water_wet_reference = [0, 0.54, 0.87, 1.46, 2.5, 3.5]
rf_water_wet_reference = [0, 0.56, 0.58, 0.60, 0.62, 0.625]

Wid_mixed_wet_reference = [0, 0.26, 0.66, 1.11, 1.89, 2.5, 3.5]
rf_mixed_wet_reference = [0, 0.39, 0.45, 0.47, 0.50, 0.51, 0.515]

Wid_plot = []
for i in Wid:
    if i < 4:
        Wid_plot.append(i)
Wid_plot.insert(0, 0)

rf_plot = []
for i in range(len(Wid_plot)):
    rf_plot.append(rf[i])
    
# Adding production data 
oiip = dataset.iloc[20, 13]
bo = dataset.iloc[21, 13]
pv = dataset.iloc[22, 13]

Npd_prod = dataset.iloc[25:, 12]
Wid_prod = dataset.iloc[25:, 13]

RF_production_ = Npd_prod*bo/oiip
Pvi_production_ = Wid_prod/pv

    
"-------------------------------------------------------------------------------------------------"   
"-------------------------------------------------------------------------------------------------"   
"-------------------------------------------------------------------------------------------------"   

"""                                   Visualising the results                                   """
                                                                   
plt.close('all')

#Water injected in PV vs RF
figure_9 = plt.figure()
plt.plot(Wid_oil_wet_reference, rf_oil_wet_reference, color = 'green', lw = 0.5, label = 'Oil-Wet ref')
plt.plot(Wid_water_wet_reference, rf_water_wet_reference, color = 'blue', lw = 0.5, label = 'Water-Wet ref')
plt.plot(Wid_mixed_wet_reference, rf_mixed_wet_reference,color = 'orange', lw = 0.5, label = 'Mixed-Wet ref')
plt.plot(Wid_plot, rf_plot, color = 'midnightblue', lw = 1.5, ls = '-', label = 'Simulated data')
plt.plot(Wid_plot, rf_plot, color = 'cornsilk', ls = '--', lw = 1.5)
plt.plot(Pvi_production_, RF_production_, color = 'red', lw = 1.5, ls = '-', label = 'Observed data')
plt.plot(Pvi_production_, RF_production_, color = 'cornsilk', ls = '--', lw = 1.5, label = '')
plt.xlabel('Water injected PV')
plt.ylabel('Recovery factor')
plt.title('Water injected PV vs RF')
plt.legend(shadow = True, loc='lower right', facecolor = 'lemonchiffon', fancybox = True, frameon=True)
plt.show()
            
"------------------------------------------------------------------------------------------------"  

krw_WW_upper =  np.load('local_data/krw_WW_upper.npy')
krow_WW_upper =  np.load('local_data/krow_WW_upper.npy')
krw_WW_lower =  np.load('local_data/krw_WW_lower.npy')
krow_WW_lower =  np.load('local_data/krow_WW_lower.npy')
Sw_OW = np.load('local_data/Sw_OW.npy')
Sw_WW = np.load('local_data/Sw_WW.npy')
krw_OW_upper =  np.load('local_data/krw_OW_upper.npy')
krow_OW_upper =  np.load('local_data/krow_OW_upper.npy')
krw_OW_lower =  np.load('local_data/krw_OW_lower.npy')
krow_OW_lower =  np.load('local_data/krow_OW_lower.npy')

# Relative perms
fig = plt.figure(figsize = (11,6))
ax = fig.add_subplot(223)
ax.plot(Sg, krg, 'b', label = 'krg' )
ax.plot(Sg, krog, 'r', label = 'krog')
ax.set_xlabel('Gas saturation')
ax.set_ylabel('kr')
ax.set_ylim(0,1)
ax.set_xlim(0,1)
ax.legend(shadow = True, loc='upper center', facecolor = 'lemonchiffon', fancybox = True, frameon=True)

ax = fig.add_subplot(221)
ax.plot(Sw, krw, 'b', label = 'krw')
ax.plot(Sw, krow, 'g', label = 'krow')
ax.set_xlabel('Water saturation')
ax.set_ylabel('kr')
ax.set_ylim(0,1)
ax.set_xlim(0,1)
ax.legend(shadow = True, loc='upper center', facecolor = 'lemonchiffon', fancybox = True, frameon=True)

ax = fig.add_subplot(122)
ax.plot(Sw_WW, krw_WW_upper, color = 'blue', lw = 0.2)
ax.plot(Sw_WW, krow_WW_upper, color = 'blue', lw = 0.2)
ax.plot(Sw_WW, krw_WW_lower, color = 'blue', lw = 0.2)
ax.plot(Sw_WW, krow_WW_lower, color = 'blue', lw = 0.2)
plt.fill_between(Sw_WW, krw_WW_upper, krw_WW_lower, alpha = 0.1, color = 'blue')
plt.fill_between(Sw_WW, krow_WW_upper, krow_WW_lower, alpha = 0.1, color = 'blue')
ax.plot(Sw_OW, krw_OW_upper, color = 'green', lw = 0.2)
ax.plot(Sw_OW, krow_OW_upper, color = 'green', lw = 0.2)
ax.plot(Sw_OW, krw_OW_lower, color = 'green', lw = 0.2)
ax.plot(Sw_OW, krow_OW_lower, color = 'green', lw = 0.2)
plt.fill_between(Sw_OW, krw_OW_upper, krw_OW_lower, alpha = 0.2, color = 'green')
plt.fill_between(Sw_OW, krow_OW_upper, krow_OW_lower, alpha = 0.2, color = 'green')
ax.plot(Sw, krw, color = 'deepskyblue', label = 'krow', linewidth = 2.5)
ax.plot(Sw, krow, color = 'lawngreen', label = 'krw', linewidth = 2.5)
ax.set_xlabel('Water saturation')
ax.set_ylabel('kr')
ax.set_ylim(0,1)
ax.set_xlim(0,1)
ax.legend(shadow = True, loc='upper center', facecolor = 'lemonchiffon', fancybox = True, frameon=True)
plt.tight_layout()
plt.show()



"-------------------------------------------------------------------------------------------------"  

# Fractional flow as a function of viscosities 
figure_3 = plt.figure()
for i, x in enumerate(fw_1):
    plt.plot(Sw_fw, x, label = str(viscosities_name[i]) + str('cp'))
    plt.xlabel('Water saturation')
    plt.ylabel('fw')
    plt.legend(shadow = True, facecolor = 'lemonchiffon', fancybox = True, frameon=True)
    plt.title('fw as a function of oil viscosities')

figure_3.tight_layout()

"-------------------------------------------------------------------------------------------------"  

# Water injected vs. Oil produced in PV
# Fractional flow 
# Water injected vs average water saturation
# Oil produced in PV vs fractional flow
plt.figure(figsize=(8,6))

plt.subplot(321)
plt.semilogx(Wid, Npd)
plt.xlabel('Water injected (PV)')
plt.ylabel('Oil produced (PV)') 
plt.ylim(0,1)
plt.grid(True)
plt.tight_layout()

plt.subplot(122)
plt.plot(x_tangent, y_tangent, color = 'black')
plt.plot(Sw_fw_plot, fw_plot, color = 'blue')
plt.scatter((Sw_at_dfw_dSw[0]), (fw_at_dfw_dSw[0]), color = 'red') 
plt.xlabel('Water saturation')
plt.ylabel('fw')
plt.ylim(0, 1)
plt.tight_layout()
plt.grid(True)

plt.subplot(323)
plt.semilogx(Wid, Sw_average)
plt.xlabel('Water injected (PV) after bt')
plt.ylabel('Average water saturation')
plt.ylim(0, 1)
plt.tight_layout()
plt.grid(True)

plt.subplot(325)
plt.plot(Npd, wc)
plt.xlabel('Oil produced (PV) after bt')
plt.ylabel('Water cut')
plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.4,
                    wspace=0.35)
a = float(Npd[-1]+0.1)
plt.xlim(0, a)
plt.grid(True)
plt.tight_layout()
plt.show()

"-------------------------------------------------------------------------------------------------"
"-------------------------------------------------------------------------------------------------"
"-------------------------------------------------------------------------------------------------"
"-------------------------------------------------------------------------------------------------"

''' REFERENCES '''


"""
https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.plot.html#matplotlib.pyplot.plot
https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D.set_transform
https://matplotlib.org/3.1.1/tutorials/introductory/pyplot.html
https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/legend_demo.html#sphx-glr-gallery-text-labels-and-annotations-legend-demo-py
https://stackoverflow.com/questions/31540258/how-to-get-alternating-colours-in-dashed-line-using-matplotlib
"""


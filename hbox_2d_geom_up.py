### getting hboxes in 2D ###

### for one cell first: setup is large enough grid, barrier that crosses it at any angle
### task: get one hbox normal to barrier on one side

# import setup:
import numpy as np
import matplotlib.pyplot as plt
import math

# physical grid setup: (also input for making function later)
xe = np.linspace(0,1,num=11)  # x grid edges
dx = xe[1]-xe[0]
xc = xe + dx/2
xc = xc[:-1]  # x grid centers
ye = np.linspace(0,1,num=11) # y grid edges
dy = ye[1] - ye[0]
yc = ye + dy/2
yc = yc[:-1] # y grid centers

# barrier parameters:
wall_height = 2 # above sea level
# indices and coordinates for start of barrier and end of barrier (x_i,y_j)
i_0 = 2 # input for making function later
j_0 = 6 # input
i_e = 7 # input
j_e = 1 # input
i_0_p = i_0
i_0_t = i_0
j_0_p = j_0
j_0_t = j_0
x_0 = xe[i_0]+0.1*dx # input for making function later
y_0 = ye[j_0]+0.3*dy # input
x_e = xe[i_e]+0.2*dx # input
y_e = ye[j_e]+0.5*dy # input


# get length of barrier:
L = np.sqrt((x_e-x_0)**2 + (y_e-y_0)**2)
theta = np.arctan((y_e-y_0)/(x_e-x_0))
m = math.ceil(L/dx) # total number of single hboxes on one side of barrier

# grid on the barrier:
# the edges (i.e. top vertices of hboxes)
we = np.zeros((2,m+1))
wc = np.zeros((2,m))
for i in range(m):
    we[0,i] = x_0 + i*np.cos(theta)*dx
    we[1,i] = y_0 + i*np.sin(theta)*dy
we[0,-1] = x_e
we[1,-1] = y_e
for i in range(m):
    wc[0,i] = 0.5*(we[0,i+1]+we[0,i])
    wc[1,i] =0.5*(we[1,i+1]+we[1,i])


# vertices of hboxes, loop by number of hboxes
ver_hbox = np.zeros((2,2*(m+1))) #fewer than 4*m because one is joined to next one
count = range(2*(m+1))
even = [num for num in count if num % 2 == 0]
odd = [num for num in count if num % 2 != 0]
ver_hbox[:,np.asarray(even)] = we
odd = np.asarray(odd)
even = np.asarray(even)
for i in range(m+1):

    ver_hbox[0,2*i+1] = we[0,i]-dx*np.sin(theta)
    ver_hbox[1,2*i+1] = we[1,i]+dy*np.cos(theta)

d1=[(ver_hbox[0,1]-ver_hbox[0,0]),(ver_hbox[1,1]-ver_hbox[1,0])]
d2 = [(ver_hbox[0,2]-ver_hbox[0,0]),(ver_hbox[1,2]-ver_hbox[1,0])]
print(np.dot(d1,d2))

# upper first edge:
we_1 = ver_hbox[:,odd]
wc_1 = np.zeros((2,m))
for i in range(m):
    wc_1[0,i] = 0.5*(we_1[0,i+1]+we_1[0,i])
    wc_1[1,i] =0.5*(we_1[1,i+1]+we_1[1,i])
# hbox centers:
hbox_cen = np.zeros((2,m))
for i in range(m):
    hbox_cen[0,i] = 0.5*(wc[0,i]+wc_1[0,i])
    hbox_cen[1,i] = 0.5*(wc_1[1,i]+wc[1,i])

# find intersection of hbox edges with physical grid edges
vertical = False
if x_e-x_0 != 0:
    slope_bar = (y_e-y_0)/(x_e-x_0)
else:
    vertical = True

# index and coordinate for start of upper edge of hbox one h-dist away
x_1= ver_hbox[0,1]
y_1 = ver_hbox[1,1]
i_1 = i_0 + int(np.sign(slope_bar)*int(xe[i_0]>x_1)) + int(np.sign(slope_bar)*int(xe[i_0+1]<x_1))
j_1 = j_0  - int(ye[j_0]>y_1) + int(ye[j_0+1]<y_1)
print("HEREE",i_1,j_1)

if slope_bar < 0:
    ## for normal edges
    i_chk = np.asarray([0,1])
    j_chk = np.asarray([0,1])

    intersect=[]
    for k in range(m+1):
        xrange = sorted([ver_hbox[0,2*k] , ver_hbox[0,2*k+1]], key = lambda x:float(x))
        yrange = sorted([ver_hbox[1,2*k] , ver_hbox[1,2*k+1]], key = lambda x:float(x))
        #(j_0)
        i_range = i_0_t + i_chk
        j_range = j_0_t + j_chk
        if any(index >= len(xe)  for index in i_range):
            #("here")
            i_range = np.asarray([len(xe)-2, len(xe)-1])
        if any(index < 0 for index in j_range) and slope_bar < 0:
            j_range = np.asarray([0,1])
        elif any(index >= len(ye) for index in j_range) and slope_bar > 0 :
            j_range = np.asarray([len(ye)-2,len(ye)-1])
        for index in i_range:
            #(xe[index])
            if xe[index] < xrange[1] and xe[index] > xrange[0]:
                dist_x = xe[index] - xrange[0]
                dist_y = dist_x * np.tan(theta-np.pi/2)
                intersect.append([xe[index],yrange[0]+dist_y])

        for index in j_range:
            #(ye[index])
            if ye[index] < yrange[1] and ye[index] > yrange[0]:
                dist_y = ye[index] - yrange[0]
                dist_x =  dist_y /np.tan(theta-np.pi/2)
                intersect.append([xrange[0]+dist_x, ye[index]])



        if k == m:
            break
        i_0_t += int(xe[i_0_t+1]<ver_hbox[0,2*(k+1)])
        j_0_t += int(np.sign(slope_bar)*int(ye[j_0_t+1]<ver_hbox[1,2*(k+1)])) + int(np.sign(slope_bar)*int(ye[j_0_t]>ver_hbox[1,2*(k+1)]))
        print("LOOK",i_0_t,j_0_t)
    ## for parallel edges:
    for k in range(m):
        xrange = sorted([ver_hbox[0,2*k] , ver_hbox[0,2*k+2]], key = lambda x:float(x))
        xrange2 = sorted([ver_hbox[0,2*k+1] , ver_hbox[0,2*k+3]], key = lambda x:float(x))
        yrange = sorted([ver_hbox[1,2*k] , ver_hbox[1,2*k+2]], key = lambda x:float(x))
        yrange2 = sorted([ver_hbox[1,2*k+1] , ver_hbox[1,2*k+3]], key = lambda x:float(x))

        #(j_0_p)
        i_range = i_0_p + i_chk
        j_range = j_0_p + j_chk
        if any(index >= len(xe)  for index in i_range):
            #("here")
            i_range = np.asarray([len(xe)-2, len(xe)-1])
        if any(index < 0 for index in j_range) and slope_bar < 0:
            j_range = np.asarray([0,1])
        elif any(index >= len(ye) for index in j_range) and slope_bar > 0 :
            j_range = np.asarray([len(ye)-2,len(ye)-1])
        for index in i_range:
            #(xe[index])
            if xe[index] < xrange[1] and xe[index] > xrange[0]:
                dist_x = xrange[1] - xe[index]
                dist_y = dist_x * np.tan(np.pi-theta)
                intersect.append([xe[index],yrange[0]+dist_y])

        for index in j_range:
            #(ye[index])
            if ye[index] < yrange[1] and ye[index] > yrange[0]:
                dist_y = yrange[1] - ye[index]
                dist_x =  dist_y /np.tan((np.pi-theta))
                intersect.append([xrange[0]+dist_x, ye[index]])


        i_0_p += int(xe[i_0_p+1]<ver_hbox[0,2*(k+1)])
        j_0_p += int(np.sign(slope_bar)*int(ye[j_0_p+1]<ver_hbox[1,2*(k+1)])) + int(np.sign(slope_bar)*int(ye[j_0_p]>ver_hbox[1,2*(k+1)]))
        print("top",i_0_p,j_0_p)

        i_range2 = i_1 + i_chk
        j_range2 = j_1 + j_chk
        if any(index >= len(xe)  for index in i_range2):
            #("here")
            i_range2 = np.asarray([len(xe)-2, len(xe)-1])
        if any(index < 0 for index in j_range2) and slope_bar < 0:
            j_range2 = np.asarray([0,1])
        elif any(index >= len(ye) for index in j_range2) and slope_bar > 0 :
            j_range2 = np.asarray([len(ye)-2,len(ye)-1])
        for index in i_range2:
            #(xe[index])
            if xe[index] < xrange2[1] and xe[index] > xrange2[0]:
                dist_x2 = xrange2[1] - xe[index]
                dist_y2 = dist_x2 * np.tan(np.pi-theta)
                intersect.append([xe[index],yrange2[0]+dist_y2])

        for index in j_range2:
            #(ye[index])
            if ye[index] < yrange2[1] and ye[index] > yrange2[0]:
                dist_y2 = yrange2[1] - ye[index]
                dist_x2 =  dist_y2 /np.tan((np.pi-theta))
                intersect.append([xrange2[0]+dist_x2, ye[index]])
        if k == m-1:
            break
        i_1 += int(xe[i_1+1]<ver_hbox[0,2*k+3])
        j_1 += int(np.sign(slope_bar)*int(ye[j_1+1]<ver_hbox[1,2*k+3])) + int(np.sign(slope_bar)*int(ye[j_1]>ver_hbox[1,2*k+3]))
        print("bottom",i_1,j_1)
    ## phys vertices inside hbox:
    i_h = i_0 + int(np.sign(slope_bar)*int(xe[i_0]>hbox_cen[0,0])) + int(np.sign(slope_bar)*int(xe[i_0+1]<hbox_cen[0,0]))
    j_h = j_0  - int(ye[j_0]>hbox_cen[1,0]) + int(ye[j_0+1]<hbox_cen[1,0])
    print(i_h,j_h)
    for k in range(m):
        cands = np.zeros((2,4))
        cands[0,:] = np.asarray([xe[i_h],xe[i_h],xe[i_h+1],xe[i_h+1]])
        cands[1,:] = np.asarray([ye[j_h],ye[j_h+1],ye[j_h],ye[j_h+1]])
        cen_x = hbox_cen[0,k]
        cen_y = hbox_cen[1,k]
        # rhombus corners
        A = ver_hbox[:,2*k]
        B = ver_hbox[:,2*k+1]
        D = ver_hbox[:,2*k+2]
        C = ver_hbox[:,2*k+3]
        a = 0.5*np.sqrt((A[0]-C[0])**2 + (A[1]-C[1])**2)
        b = 0.5*np.sqrt((B[0]-D[0])**2 + (B[1]-D[1])**2)
        U = (C-A)/(2*a)
        V = (D-B)/(2*b)
        Q = np.asarray([cen_x,cen_y])
        for l in range(4):
            P = np.asarray([cands[0,l],cands[1,l]])
            W = P-Q
            if abs(np.dot(W,U)) /a + abs(np.dot(W,V)) / b <= 1:
                intersect.append([cands[0,l],cands[1,l]])
                print(l, cen_x,cen_y,cands[0,l],cands[1,l])
        if k==m-1:
            break
        i_h += -int(xe[i_h]>hbox_cen[0,k+1]) +int(xe[i_h+1]<hbox_cen[0,k+1])
        j_h +=int(ye[j_h+1]<hbox_cen[1,k+1])- int(ye[j_h]>hbox_cen[1,k+1])
        print(i_h,j_h)
    last_mid_x = 0.5*(ver_hbox[0,-2] + ver_hbox[0,-1])
    last_mid_y = 0.5*(ver_hbox[1,-2] + ver_hbox[1,-1])
    intersection = np.asarray(intersect)
    print(last_mid_x,intersection[-1,0])

    if intersection[-1,0] > last_mid_x:
        print("here!!!!!")
        intersect.pop()
        intersection = np.asarray(intersect)
print(slope_bar)

##################################################33
if slope_bar > 0:

    print("here")
    ###########################################
    ##                                       ##
    ##                                       ##
    ## for normal edges                      ##
    i_chk = np.asarray([-1,0,1])
    j_chk = np.asarray([-1,0,1])

    intersect=[]
    for k in range(m+1):
        xrange = sorted([ver_hbox[0,2*k] , ver_hbox[0,2*k+1]], key = lambda x:float(x))
        yrange = sorted([ver_hbox[1,2*k] , ver_hbox[1,2*k+1]], key = lambda x:float(x))

        #(j_0)
        i_range = i_0_t + i_chk
        j_range = j_0_t + j_chk

        if any(index >= len(xe)  for index in i_range):
            #("here")
            i_range = np.asarray([len(xe)-2, len(xe)-1])
        if any(index < 0 for index in j_range) and slope_bar < 0:
            j_range = np.asarray([0,1])
        elif any(index >= len(ye) for index in j_range) and slope_bar > 0 :
            j_range = np.asarray([len(ye)-2,len(ye)-1])
        for index in i_range:
            print([xrange[0],xe[index],xrange[1]])
            if xe[index] < xrange[1] and xe[index] > xrange[0]:
                print("ME",xe[index])
                dist_x = abs(xrange[0] - xe[index])
                dist_y = abs(dist_x /np.tan(theta))
                intersect.append([xe[index],yrange[1]-dist_y])

        for index in j_range:
            #(ye[index])
            if ye[index] < yrange[1] and ye[index] > yrange[0]:
                dist_y = abs(yrange[1] - ye[index])
                dist_x =  abs(dist_y/np.tan(np.pi/2-theta))
                intersect.append([xrange[0]+dist_x, ye[index]])


        if k == m:
            break
        i_0_t += int(xe[i_0_t+1]<ver_hbox[0,2*(k+1)])
        j_0_t += int(np.sign(slope_bar)*int(ye[j_0_t+1]<ver_hbox[1,2*(k+1)])) + int(np.sign(slope_bar)*int(ye[j_0_t]>ver_hbox[1,2*(k+1)]))
        print("I,J",i_0_t,j_0_t)
        print("IRANGE,JRANGE",i_range,j_range)
    ##################################
    ## for parallel edges:
    ###################################
    for k in range(m):
        xrange = sorted([ver_hbox[0,2*k] , ver_hbox[0,2*k+2]], key = lambda x:float(x))
        xrange2 = sorted([ver_hbox[0,2*k+1] , ver_hbox[0,2*k+3]], key = lambda x:float(x))
        yrange = sorted([ver_hbox[1,2*k] , ver_hbox[1,2*k+2]], key = lambda x:float(x))
        yrange2 = sorted([ver_hbox[1,2*k+1] , ver_hbox[1,2*k+3]], key = lambda x:float(x))

        #(j_0_p)
        i_range = i_0_p + i_chk
        j_range = j_0_p + j_chk
        if any(index >= len(xe)  for index in i_range):
            #("here")
            i_range = np.asarray([len(xe)-2, len(xe)-1])
        if any(index < 0 for index in j_range) and slope_bar < 0:
            j_range = np.asarray([0,1])
        elif any(index >= len(ye) for index in j_range) and slope_bar > 0 :
            j_range = np.asarray([len(ye)-2,len(ye)-1])
        for index in i_range:
            #(xe[index])
            if xe[index] < xrange[1] and xe[index] > xrange[0]:
                dist_x = abs(xrange[0] - xe[index])
                dist_y = dist_x * np.tan(theta)
                intersect.append([xe[index],yrange[0]+dist_y])

        for index in j_range:
            #(ye[index])
            if ye[index] < yrange[1] and ye[index] > yrange[0]:
                dist_y = abs(yrange[0] - ye[index])
                dist_x =  dist_y /np.tan((theta))
                intersect.append([xrange[0]+dist_x, ye[index]])


        i_0_p += int(xe[i_0_p+1]<ver_hbox[0,2*(k+1)])
        j_0_p += int(np.sign(slope_bar)*int(ye[j_0_p+1]<ver_hbox[1,2*(k+1)])) + int(np.sign(slope_bar)*int(ye[j_0_p]>ver_hbox[1,2*(k+1)]))
        print("bottom",i_0_p,j_0_p)
        i_range2 = i_1 + i_chk
        j_range2 = j_1 + j_chk
        if any(index >= len(xe)  for index in i_range2):
            #("here")
            i_range2 = np.asarray([len(xe)-2, len(xe)-1])
        if any(index < 0 for index in j_range2) and slope_bar < 0:
            j_range2 = np.asarray([0,1])
        elif any(index >= len(ye) for index in j_range2) and slope_bar > 0 :
            j_range2 = np.asarray([len(ye)-2,len(ye)-1])
        for index in i_range2:
            #(xe[index])
            if xe[index] < xrange2[1] and xe[index] > xrange2[0]:
                dist_x2 = abs(xrange2[0] - xe[index])
                dist_y2 = dist_x2 * np.tan(theta)
                intersect.append([xe[index],yrange2[0]+dist_y2])

        for index in j_range2:
            #(ye[index])
            if ye[index] < yrange2[1] and ye[index] > yrange2[0]:
                dist_y2 = abs(yrange2[0] - ye[index])
                dist_x2 =  dist_y2 /np.tan((theta))
                intersect.append([xrange2[0]+dist_x2, ye[index]])
        if k == m-1:
            break
        i_1 += int(xe[i_1+1]<ver_hbox[0,2*k+3])
        j_1 += int(np.sign(slope_bar)*int(ye[j_1+1]<ver_hbox[1,2*k+3])) + int(np.sign(slope_bar)*int(ye[j_1]>ver_hbox[1,2*k+3]))
        print("top",i_1,j_1)

    ## phys vertices inside hbox:
    i_h = i_0 + int(np.sign(slope_bar)*int(xe[i_0]>hbox_cen[0,0])) + int(np.sign(slope_bar)*int(xe[i_0+1]<hbox_cen[0,0]))
    j_h = j_0  - int(ye[j_0]>hbox_cen[1,0]) + int(ye[j_0+1]<hbox_cen[1,0])
    print(i_h,j_h)
    for k in range(m):
        cands = np.zeros((2,4))
        cands[0,:] = np.asarray([xe[i_h],xe[i_h],xe[i_h+1],xe[i_h+1]])
        cands[1,:] = np.asarray([ye[j_h],ye[j_h+1],ye[j_h],ye[j_h+1]])
        cen_x = hbox_cen[0,k]
        cen_y = hbox_cen[1,k]
        # rhombus corners
        A = ver_hbox[:,2*k]
        B = ver_hbox[:,2*k+1]
        D = ver_hbox[:,2*k+2]
        C = ver_hbox[:,2*k+3]
        a = 0.5*np.sqrt((A[0]-C[0])**2 + (A[1]-C[1])**2)
        b = 0.5*np.sqrt((B[0]-D[0])**2 + (B[1]-D[1])**2)
        U = (C-A)/(2*a)
        V = (D-B)/(2*b)
        Q = np.asarray([cen_x,cen_y])
        for l in range(4):
            P = np.asarray([cands[0,l],cands[1,l]])
            W = P-Q
            if abs(np.dot(W,U)) /a + abs(np.dot(W,V)) / b <= 1:
                intersect.append([cands[0,l],cands[1,l]])
                print(l, cen_x,cen_y,cands[0,l],cands[1,l])
        if k==m-1:
            break
        i_h += -int(xe[i_h]>hbox_cen[0,k+1]) +int(xe[i_h+1]<hbox_cen[0,k+1])
        j_h +=int(ye[j_h+1]<hbox_cen[1,k+1])- int(ye[j_h]>hbox_cen[1,k+1])
        print(i_h,j_h)
    last_mid_x = 0.5*(ver_hbox[0,-2] + ver_hbox[0,-1])
    last_mid_y = 0.5*(ver_hbox[1,-2] + ver_hbox[1,-1])
    intersection = np.asarray(intersect)
    print(last_mid_x,intersection[-1,0])

    if intersection[-1,0] > last_mid_x:
        print("here!!!!!")
        intersect.pop()
        intersection = np.asarray(intersect)
print(intersection)








## PLOTTING ##
plt.figure()
xeg, yeg = np.meshgrid(xe,ye)
plt.plot(xeg,yeg,'k:')
plt.plot(yeg,xeg,'k:')
plt.plot(np.linspace(x_0,x_e), np.tan(theta)*np.linspace(x_0,x_e) + (y_e - np.tan(theta)*x_e),'r')
plt.plot(ver_hbox[0,np.asarray(odd)],ver_hbox[1,np.asarray(odd)],'c')
for i in range(m+1):
    plt.plot(ver_hbox[0,[2*i,2*i+1]],ver_hbox[1,[2*i,2*i+1]],'c')
for k in range(len(intersection[:,1])):
    plt.plot(intersection[k,0],intersection[k,1],'g*')
plt.plot(wc[0,:],wc[1,:],'k*')
plt.plot(wc_1[0,:],wc_1[1,:],'k*')
plt.plot(hbox_cen[0,:],hbox_cen[1,:],'k*')
plt.show()

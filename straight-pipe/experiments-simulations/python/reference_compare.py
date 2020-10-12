#! /usr/bin/python
# -*- coding: utf-8 -*-

# this file is for reference Fr data without axes and fix the bubble position
# read out data auto

from __future__ import unicode_literals
import os
import sys
import numpy as np
import math
import matplotlib 
import matplotlib.pyplot as plt
import pylab
import argparse
import string
import math
import re
import operator
from scipy import ndimage
#import pyautogui


from PIL import Image
from PIL import ImageEnhance
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['figure.dpi'] = 50
matplotlib.rcParams['figure.dpi'] = 50

CWP      = os.path.abspath("./")


case = [16,28,33,47,105,290]
case_1 = [1,2,3,4,5,6,7,8,9,10]
Eo = [6.34,9,10.54,15.25,18.67,42.4,75.44,131.84,300]
#Eo = [5.5,7.9,9.04,13.04,15.8,34.94,106.95,150,111,70]
#for k in case:
#    for m in case_1:


a = 8
b = 0
while b < 5:
    #k = int(input('Please input Mo number, choose from (16,28,33.47,105,290):'))  
    while a<10:
        i = b
        a = a + 1
        try:
            #a = int(input('Please input the case number from 1 to 8:'))
            j = case_1.index(a)

            PATH_FR = os.path.abspath('/media/lguo/Elements/phd-manuscript/test-cases/no-phase-change/taylor-bubble/2D/straight_pipe/reference-paper/{}_simulation/{}'.format(case[i],case_1[j]))
            if not os.path.exists(PATH_FR):
                os.makedirs(PATH_FR)
            PATH_FR_DATA = PATH_FR + '/taylor_2D'
            PATH_INTERFACE = PATH_FR_DATA + '/interface'
            PATH_VELOCITY =  PATH_FR_DATA + '/velocity'

            REFERENCE_SAVE = os.path.abspath('/media/lguo/Elements/phd-manuscript/thesis-manuscript/Taylor-bubble/latex/elsarticle-taylor/figure/straight/{}/{}-{}'.format(case[i],case[i],case_1[j]))
            if not os.path.exists(REFERENCE_SAVE):
                os.makedirs(REFERENCE_SAVE)


            def get_interface(time):
                fn = os.path.join(PATH_INTERFACE, 'interface-%g' % time)
                sur = []
                with open(fn, 'r') as file:
                    seg = []
                    for line in file.readlines():
                        line = line.rstrip('\n')
                        l = line.split(' ')
                        if len(l) == 2:
                            #print float(l[0]), float(l[1])
                            seg.append([float(l[0]), float(l[1])])
                        else:
                            segn = seg[:]
                            sur.append(segn)
                            seg = []
                return sur

            def get_velocity(time):
                fn = os.path.join(PATH_VELOCITY, 'vprof-%.4f' % time)
                arrs = []
                strr = ''
                with open(fn, 'r') as file:
                    for line in file:
                        if '1e+30' in line:
                            line=line.replace('1e+30','0')
                        strr += line
                with open(fn, 'w') as file:
                    file.write(strr)
                with open(fn, 'r') as file:
                    for line in file.readlines():
                        line = line.strip().strip('\n')
                        l = line.split(' ')
                        #print float(l[0]), float(l[1]),float(l[2]),float(l[3])
                        arrs.append([float(l[0]), float(l[1]),float(l[2]),float(l[3])])
                return arrs

            def get_bubble_head(sur):
                maxx = sur[0][0][0]
                for seg in sur:
                    x1 = seg[0][0]
                    y1 = seg[0][1]
                    x2 = seg[1][0]
                    y2 = seg[1][1]
                    if x1 > maxx:
                        maxx = x1
                    if x2 > maxx:
                        maxx = x2
                return maxx

            def get_bubble_tail(sur):
                minn = sur[0][0][0]
                for seg in sur:
                    x1 = seg[0][0]
                    y1 = seg[0][1]
                    x2 = seg[1][0]
                    y2 = seg[1][1]
                    if x1 < minn:
                        minn = x1
                    if x2 < minn:
                        minn = x2
                return minn

            def out_data(filename):
                fn = os.path.join(PATH_FR_DATA , filename)
                seg = []
                with open(fn, 'r') as file:
                    for line in file.readlines()[0:2500]:
                        line = line.rstrip('\n')
                        l = line.split(' ')
                        if len(l) == 26:
                            #print float(l[0]), float(l[1])
                            seg.append([float(l[0]), float(l[5]), float(l[6])])
                        else:
                            seg = []
                return seg

            def plot_sur(plt, sur, sig, color):
                
                for seg in sur:
                    x1 = seg[0][0]
                    x2 = seg[1][0]
                    y1 = seg[0][1]
                    y2 = seg[1][1]
                    plt.plot([y1, y2], [x1 - sig, x2 - sig], color,lw = 1)
                    res, = plt.plot([-y1, -y2], [x1- sig, x2 - sig], color, lw = 1)
                return res


            def plot_velocity(plt, sur): 
                soa =np.array(sur) 
                X,Y,U,V = zip(*soa)
                M  = np.hypot(U, V)
                ax = plt.gca()
                YN = np.array([-val for val in Y])
                ax.quiver(YN,X,V,U,units = 'xy', angles='xy',scale_units='xy',  scale = 10,width = 0.005,color = 'k')


            def plot_legend(plt,time, xmin, ymin, xmax, ymax):

                str_undim = 'log(Mo) = -4.63\nEo =40.55\nNf = 230.9\nRatio=910'
                plt.axes().text(xmin - 2.5, (ymin+0.5), 
                    str_undim,
                    bbox={'facecolor':'white'}, 
                    fontname='monospace',
                    )
                str_time = r'$t^*$ = %.2f' % (time) 
                plt.axes().text(xmin - 2, (ymin-0.1), 
                    str_time,bbox={'facecolor':'white'}, 
                    fontname='monospace',
                    color='r'
                    )

            def plot_streamline(plt,sur,a,time,ex,tail,head):
                from scipy.interpolate  import griddata
                ax = plt.gca()
                soa =np.array(sur) 
                x,y,u,v = zip(*soa)
                nx, ny = 2000, 2000
                pts = np.vstack((x, y)).T
                vals = np.vstack((u, v)).T
                xi = np.linspace(0,16,2000)
                yi = np.linspace(0,0.5,2000)
                ipts = np.vstack(a.ravel() for a in np.meshgrid(yi, xi)[::-1]).T
                ivals = griddata(pts, vals, ipts, method='cubic')
                ui, vi = ivals.T
                ui.shape = vi.shape = (ny, nx)
                time_1 = time*100 + 1
                b = int(time_1)
                plt.axvline(0,ymin = -3, ymax = 3, color='b',linestyle = '-.',linewidth = 0.5)
                #stream_points_1 = np.array(zip(np.arange(0,0.5,.05), np.arange(tail-2,head+0.5,.1)))
                stream_points = np.array(zip(np.arange(0,0.5,.05), np.arange(head-2,head,.1)))
                stream_points_1 = np.array(zip(np.arange(0,0.5,.05), np.arange(head-3.5,head-3.3,.1)))
                ax.streamplot(yi, xi, vi-a[b][2], ui-a[b][1], color = 'k',start_points=stream_points, arrowsize = 0.5,linewidth = 0.5,density =50)
                ax.streamplot(yi, xi, vi-a[b][2], ui-a[b][1], color = 'k',start_points=stream_points_1,arrowsize = 0.5,linewidth = 0.5,density =50)

            def simulation():

                for time in np.arange(7,8,1):

                    plt.figure(figsize=(3, 5))

                    """
                    Set labels
                    """
                    #plt.xlabel(r'$r/D$')
                    #plt.ylabel(r'$(z_h - z)/D$')
                    plt.axis('off')
                    """
                    Set range
                    """
                    #plt.xscale('log')
                    #plt.yscale('log')
                    # plt.ylim([0,0.5])

                    sur1 = get_interface(time)
                    sur2 = get_velocity(time)
                    plot_velocity(plt,sur2)
                    a = out_data('out')

                    ex= [1.72,1.45,1.33,1.24,1.12,1,0.9,0.81,0.69] # expansion ratio
                    
                    xh1 = get_bubble_head(sur1)
                    xh11 = get_bubble_tail(sur1)
                    print "head tail:", xh1, xh11
                    plot_streamline(plt,sur2,a,time,ex,xh11,xh1)

                    with open('data', 'w') as f:  
                        f.write("head-%.1f :%f\t" %(time, xh1))
                        f.write("tail-%.1f :%f\n" %(time, xh11))
                    le1 = plot_sur(plt, sur1, 0, "r")
                    
                    ex_ratio = ex[5]       
                    if ex_ratio in ex[:5]:
                        expansion = ex_ratio + 0.08
                        plt.xlim(-0.5*expansion,0.5*expansion)
                    else:
                        expansion = 1 + 0.08
                        plt.xlim(-0.5*expansion,0.5*expansion)

                    y_min = xh1 - 3.5
                    y_max = xh1 + 0.5
                    plt.ylim(y_min,y_max)
                    y1 = math.fabs(y_min)/(y_max - y_min)
                    x1 = math.fabs(0.5*expansion+0.5)/expansion
                    x2 = math.fabs(0.5+expansion)/expansion - 0.08
                    x3 = 0
                    x4 = (0.5 - 0.08)/expansion
                
                    plt.axvline(x = 0.5,ymin = 0, ymax = y1, color='k')
                    plt.axvline(x = -0.5,ymin = 0, ymax = y1, color='k')
                    plt.axvline(x = 0.5*ex_ratio,ymin = y1, ymax = 1, color='k')
                    plt.axvline(x = -0.5*ex_ratio,ymin = y1, ymax = 1, color='k')
                    plt.axhline(y = 0,xmin=0.08*0.5/expansion,xmax=(expansion*0.5-0.5)/expansion,color='k')
                    plt.axhline(y = 0,xmin=(0.5+0.5*expansion)/expansion,xmax=1-0.08*0.5/expansion,color='k')
                    plt.axhspan(xmin = x1,xmax = x2,ymin = y_min,ymax = 0,facecolor = 'w',alpha = 1)
                    plt.axhspan(xmin = x3,xmax = x4,ymin = y_min,ymax = 0,facecolor = 'w',alpha = 1)
                    plt.title('$Eo=$'+str(Eo[j]))
                    version = matplotlib.__version__

                    plt.axes().set_aspect('equal')
                    plt.tight_layout()

                    name =  "{}-{}".format(case[i],case_1[j]) + "-%g"% (time)
                    fn = os.path.join(REFERENCE_SAVE, (name + '.pdf'))
                    pdf = PdfPages(fn.format(version))
                    #plt.savefig(fn, format='png', bbox_inches='tight', pad_inches = 0,dpi = 300)
                    plt.savefig(fn, format='pdf', bbox_inches='tight')
                    plt.close() 


            def main():
                simulation()


            if __name__ == '__main__':
                print("-----start-----")
                main()    
                print("-----finished!!!-----")
        except IOError or ValueError:
            pass
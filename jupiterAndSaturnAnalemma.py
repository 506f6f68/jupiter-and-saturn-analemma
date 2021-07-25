import numpy as np
import matplotlib.pyplot as plt
import matplotlib; matplotlib.use('agg')
from mpl_toolkits.mplot3d import Axes3D
import imageio
import csv
import astropy.time
import astropy.time.formats
import astropy.units as u
from astropy.coordinates import solar_system_ephemeris,SkyCoord,EarthLocation,AltAz 
import os

class analemaa: 
    def __init__ (self,obs_loc, start, end, time,timedelta):
        self.obs_loc = obs_loc
        self.start = start
        self.end = end
        self.time = time
        self.timedelta = timedelta

        self.day_past = int((astropy.time.Time(self.end)-astropy.time.Time(self.start)).value)+1 # day past since start
        self.obs_time = astropy.time.Time(f'{self.start}T{self.time}')+ np.arange(self.day_past)*self.timedelta # list of observation time
        index = str(self.obs_time[0]).find('T')
        self.date_label = [str(date)[index-5:index] for date in self.obs_time]

        #Jupiter
        jup_loc = astropy.coordinates.get_body("jupiter", self.obs_time).transform_to(AltAz(obstime=self.obs_time,location=self.obs_loc)) # get jupiter position 
        self.az_j,self.alt_j = jup_loc.az.value,jup_loc.alt.value # get jupiter altitude and azimuth 
        #Saturn 
        sat_loc = astropy.coordinates.get_body("saturn", self.obs_time).transform_to(AltAz(obstime=self.obs_time,location=self.obs_loc)) # get saturn position 
        self.az_s,self.alt_s = sat_loc.az.value,sat_loc.alt.value # get saturn altitude and azimuth 
        
        self.path = os.path.dirname(os.path.abspath(__file__))

    def writeCSV (self):
        data = open(f'{self.path}/analemaa csv.csv', 'w')
        writer = csv.writer(data)
        writer.writerow(["Date", "Jupiter Altitude", "Jupiter Azimuth", "Saturn Altitude", "Saturn Azimuth"])
        for i in range (self.day_past): 
            writer.writerow([self.date_label[i],self.alt_j[i],self.az_j[i], self.alt_s[i], self.az_s[i]])

    def plotAnalemma (self):
        plt.style.use('dark_background')
        plt.scatter(self.az_j,self.alt_j,s=3, color = "red") #Plot Jupiter
        plt.scatter(self.az_s,self.alt_s,s=3, color = "yellow") #Plot Satrun 
        selectedDates = [0, round((self.day_past)/2) ,self.day_past-1] #Plot selected dates 
        for t in selectedDates:
            plt.text(self.az_j[t]+5,self.alt_j[t]+0.1,self.date_label[t],fontname='AppleGothic',ha='center', size= 8, color= "red") #Plot Jupiter Dates
            plt.text(self.az_s[t]+5,self.alt_s[t]+0.1,self.date_label[t],fontname='AppleGothic',ha='center', size= 8, color="yellow") #Plot Saturn Dates
        plt.title(f'{self.start}~{self.end} {self.time}',fontname='AppleGothic')#Plot title 
        plt.xlabel('Azimuth',fontname='AppleGothic') #Plot x label 
        plt.ylabel('Altitude',fontname='AppleGothic') #Plot y label 
        plt.savefig(f'{self.path}/ analemma.png', dpi=300, bbox_inches='tight')
    
    def analemmaAnimation (self):
        plt.style.use('dark_background')
        gif = []
        for i in range(self.day_past):
            t = i
            s = self.date_label[i]
            fig = plt.figure(figsize=[24,11])
            ax = plt.axes(xlim=[0,360],ylim=[-65,55])

            #Jupiter
            ax.scatter(self.az_j[:t],self.alt_j[:t],s=10,color = 'red')
            ax.scatter(self.az_j[t],self.alt_j[t],c='w',s=100,marker='*',edgecolor='r')

            #Saturn
            ax.scatter(self.az_s[:t],self.alt_s[:t],s=10,color ='yellow')
            ax.scatter(self.az_s[t],self.alt_s[t],c='w',s=100,marker='*',edgecolor='y')

            ax.text((self.az_j[i]+self.az_s[i])/2,(self.alt_j[i]+self.alt_s[i])/2+5,s,fontname='AppleGothic',ha='center', size= 12, color= "white")
            
            ax.set_title(f'{self.start}~{self.end} {self.time}',fontname='AppleGothic', size = 30) 
            ax.set_xlabel(u'Azimuth',fontname='AppleGothic', size = 25)
            ax.set_ylabel(u'Altitude',fontname='AppleGothic', size = 25)
            
            fig.canvas.draw()
            gif.append(np.array(fig.canvas.renderer._renderer)[65:-65,220:-190])
            
            plt.close()
        imageio.mimsave(f'{self.path}/ analemaa animation.gif',gif,fps=24)
    
    def globeAnimation (self):
        plt.style.use('default')

        #Jupiter
        x_j = -np.cos(np.radians(self.az_j))/np.tan(np.radians(self.alt_j))
        y_j = -np.sin(np.radians(self.az_j))/np.tan(np.radians(self.alt_j))
        sx_j,sy_j,sz_j = SkyCoord(ra=self.az_j,dec=self.alt_j,unit='deg').cartesian.xyz.value*2
        sz_j += 1

        #Saturn
        x_s = -np.cos(np.radians(self.az_s))/np.tan(np.radians(self.alt_s))
        y_s = -np.sin(np.radians(self.az_s))/np.tan(np.radians(self.alt_s))
        sx_s,sy_s,sz_s = SkyCoord(ra=self.az_s,dec=self.alt_s,unit='deg').cartesian.xyz.value*2
        sz_s += 1

        gif = []
        for i in range(self.day_past):
            t = i
            fig = plt.figure(figsize=[17,17])
            ax = plt.axes([0,0,1,1],projection='3d',xlim=[-1.5,1.5],ylim=[-1.5,1.5],zlim=[0,3])
            ax.invert_yaxis()

            ax.scatter(sx_j[:t],sy_j[:t],sz_j[:t],c=np.arange(t),s=10,cmap='Reds')
            ax.scatter(sx_j[t],sy_j[t],sz_j[t],c='w',s=100,marker='*',edgecolor='r')

            ax.scatter(sx_s[:t],sy_s[:t],sz_s[:t],c=np.arange(t),s=10,cmap='Wistia')
            ax.scatter(sx_s[t],sy_s[t],sz_s[t],c='w',s=100,marker='*',edgecolor='y')

            ax.text((sx_j[t]+sx_s[t])/2,(sy_j[t]+sy_s[t])/2-0.3,(sz_j[t]+sz_s[t])/2,self.date_label[t],fontname='AppleGothic')

            plt.plot([0,0],[0,0],[0,1],'k',lw=2)
            ax.view_init(15,-10)
            ax.set_yticklabels(['' for _ in ax.get_yticks()])
            fig.canvas.draw()
            gif.append(np.array(fig.canvas.renderer._renderer)[65:-65,220:-190])
            plt.close()
        
        imageio.mimsave(f'{self.path}/ globe animation.gif',gif,fps=24)
    

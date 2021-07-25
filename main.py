from jupiterAndSaturnAnalemma import*

path = os.path.dirname(os.path.abspath(__file__))

#UserInput
solar_system_ephemeris.set(f"{path}/de441_part-1.bsp") #select ephemeris (url or directory)
obs_loc = EarthLocation(lat='31.7054 deg',lon='35.2024 deg') # observe location 
start = '-00006-01-01' #observation start date
end = '-00006-12-31' #observation end date 
time = '15:00:00.000' #observation time
timedelta = 1*u.day #observation step size

#Available Function
analemaa = analemaa (obs_loc,start, end, time,timedelta)
analemaa.plotAnalemma ()
analemaa.writeCSV ()
analemaa.analemmaAnimation ()
analemaa.globeAnimation ()






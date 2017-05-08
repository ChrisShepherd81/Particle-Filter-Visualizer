'''
Created on 05.05.2017

@author: christian@inf-schaefer.de
'''

import sys, getopt
import numpy as np
import matplotlib
import math
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.animation as animation

def plotData(dataFolder, particlesFile, numberOfParticles, outputFile ):
    
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'))

    noP = int(numberOfParticles)
    #Load data
    print("Loading data ...")
    groundTruth = np.loadtxt(dataFolder + "gt_data.txt", usecols=(0, 1,2))
    mapData = np.loadtxt(dataFolder + "map_data.txt", usecols=(0, 1))
    particles = np.loadtxt(particlesFile)
    
    fig = plt.figure(figsize=(10,10))
    
    lbl_gt = mlines.Line2D([], [], color='g', label='ground truth')
    plt.legend(handles=[lbl_gt])

    #Global map
    ax = plt.subplot(211)
    plt.title('Global map', loc='left')
    fig.gca().set_ylim([-110,50])
    timeStamp_text = ax.text(-50,40,"Timestamp: 0")
    plt.plot(groundTruth[:,0], groundTruth[:,1], 'g', label='ground truth')
    plt.plot(mapData[:,0], mapData[:,1],'cD', markerfacecolor='None', label="map data")
    glb_obs, = plt.plot([],[], 'k+', markersize=8, linewidth=1.5, label="observation")
    glb_pat, = plt.plot(particles[:noP,0], particles[:noP,1], 'k.', label="particel" )
   
    #Particles view
    ax_particles = plt.subplot(223)
    plt.title('Particles', loc='left')
    ax_particles.set_xlim(-0.6,0.6)
    ax_particles.set_ylim(-0.6,0.6)
    plt.plot(0,0, 'go', markerfacecolor='None', markersize=10)
    partAni, = plt.plot([],[], 'k.', label="particle")
    part_glb, = plt.plot([],[], 'g')
    
    #Observations view
    ax_obs_view = plt.subplot(224)
    plt.title('Observations', loc='left')
    ax_obs_view.set_xlim(-60,60)
    ax_obs_view.set_ylim(-60,60)
    ax_obs_view.add_artist(plt.Circle((0, 0), 50.0, color='k', fill=False))
    patch = plt.Arrow(0, 0, 20*math.cos(0), 20*math.sin(0), linewidth=1.5, color='k')
    ax_obs_view.add_patch(patch)
    obs_view, = plt.plot([], [], 'k+', label="observations")
    obs_view_gt, = plt.plot([], [], 'cD', markerfacecolor='None', markersize=10, label="map data")
    plt.plot(0,0, 'ko', markersize=10)
    
    #legend
    ax.legend(loc=4)
    
    plt.tight_layout()
    
    with writer.saving(fig, outputFile + ".mp4", dpi=100):
        max_range = 2444
        for i in range(max_range):
            print("Timestamp: " + str(i+1))
            timeStamp_text.set_text("Timestamp: " + str(i+1))
            
            #Load observation
            observations = np.loadtxt(dataFolder + "observation/observations_" + str(i+1).zfill(6) +".txt")
            
            #Next particles for observation
            pat = particles[i*noP:(i+1)*noP,:]
            
            #Estimate car position and yaw angle
            car_pos_x = np.average(pat[:,0])
            car_pos_y = np.average(pat[:,1])
            theta = np.average(pat[:,2] )
             
            #Draw estimated theta
            ax_obs_view.patches.remove(patch)
            patch = plt.Arrow(0, 0, 20*math.cos(theta), 20*math.sin(theta), width=1.5, color='k')
            ax_obs_view.add_patch(patch)
            #Draw map markers in observation view                          
            obs_view_gt.set_data(mapData[:,0] - car_pos_x , mapData[:,1]- car_pos_y)
            
            #Draw observations in observations view
            rotMatrix = np.array([[np.cos(theta), -np.sin(theta)], \
                                  [np.sin(theta), np.cos(theta)]])  
            rot_observations = np.dot(observations, rotMatrix.T);  
            obs_view.set_data(rot_observations[:,0], rot_observations[:,1])
            
            #Draw observations in global map view
            rotMatrix = np.array([[np.cos(theta), -np.sin(theta)], \
                                   [np.sin(theta),  np.cos(theta)]])
            observations = np.dot(observations, rotMatrix.T);
            glb_obs.set_data(observations[:,0] + car_pos_x, observations[:,1] + car_pos_y)
           
            #Draw particles in global view
            glb_pat.set_data(pat[:,0], pat[:,1])
            
            #Draw particles in particle view
            gt_pos = [groundTruth[i,0], groundTruth[i,1]]
            patX = pat[:,0] - gt_pos[0]
            patY = pat[:,1] - gt_pos[1]
            partAni.set_data(patX, patY) 
            
            #Update particle view axis
            xmin, xmax = ax_particles.get_xlim()
            ymin, ymax = ax_particles.get_ylim()
            max_val = max(max(abs(patX)),max(abs(patY)), abs(xmin), abs(ymin), ymax, xmax)
            ax_particles.set_xlim(-max_val, max_val)
            ax_particles.set_ylim(-max_val, max_val)
           
            #Draw ground truth in particle view
            gt_fragment = np.array([[0,0]])
            
            idx = 0 if i < 3 else i - 3
            idx = idx if idx <= max_range - 6 else max_range - 6
            
            for j in range(0,6):
                gt_fragment = np.append(gt_fragment, [[ groundTruth[idx+j,0] - gt_pos[0], groundTruth[idx+j,1] - gt_pos[1] ]] , axis=0)
            part_glb.set_data(gt_fragment[1:,0],gt_fragment[1:,1])
                
            #Update axis
            xmin, xmax = ax_particles.get_xlim()
            ymin, ymax = ax_particles.get_ylim()
            max_val = max(max(abs(patX)),max(abs(patY)), abs(xmin), abs(ymin), ymax, xmax)
            ax_particles.set_xlim(-max_val, max_val)
            ax_particles.set_ylim(-max_val, max_val)
            writer.grab_frame()

def main(argv):
    dataFolder = ''
    particlesFile = ''
    numberOfParticles = 0
    outputFile = "animation"
    
    usage = 'ParticleFilterVisualizer.py -f <dataFolder> -p <particleFile> -n <numberOfParticles> [-o <outputFileName>]'
    
    try:
        opts, args = getopt.getopt(argv,"hf:p:n:o:",["dataFolder=","particleFile=","numberParticles=" ,"outPutFile="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
        
    if len(opts) is 0:
        print(usage)
        sys.exit()
        
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-f", "--dataFolder"):
            dataFolder = arg
        elif opt in ("-p", "--particleFile"):
            particlesFile = arg
        elif opt in ("-n", "--numberOfParticles"):
            numberOfParticles = arg
        elif opt in ("-o", "--outputFile"):
            outputFile = arg
        else:
            print(usage)
            sys.exit(2)
            
    plotData(dataFolder, particlesFile, numberOfParticles, outputFile)

if __name__ == "__main__":
    main(sys.argv[1:])
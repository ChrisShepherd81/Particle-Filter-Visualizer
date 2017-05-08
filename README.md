# Particle-Filter-Visualizer
Visualizer for the Particle-Filter

[screenshot](Screenshot.png)

## Dependencies
* python-matplotlib https://matplotlib.org/users/installing.html
* NumPy https://docs.scipy.org/doc/numpy-1.10.1/user/install.html
* FFMPEG https://ffmpeg.org/

## Running the script
Run the following commands from the command line:
```
python ParticleFilterVisualizer.py -f <dataFolder> -p <particleFile> -n <numberOfParticles> [-o <outputFileName>]

e.g.: python ParticleFilterVisualizer.py -f /data/ -p particles.txt -n 100 -o myFile
```
# mesh part
import numpy as np

Bat_to_Landscape = 8.5 # ditance to the landscape
Azimuth = 0 # head direction (depend on bats)
elevation_angle = 90 # target angle to the center point = 180 - elevation_angle
mikes_dist_coef = 0.0035 # pip kuli
# Azimuth = 0

# echo part
min_freq = 35e3
max_freq = 90e3 # for Pipistrellus kuhlii
step_size = 50  # frequency step size
freqvec = np.arange(min_freq, max_freq+1, step_size) # Emitter frequencies 25 kHz to 80 kHz
# print(freqvec)

Fs= 4096*step_size*2  # sample rate


if Fs < 3* max_freq:
	print('Warning, Raise step size')

Em_level = 1 # use 1 Pa as emission
# freq = 10 # has no function except for bemfa code
max_gain_of_ear_db = 21
mic_angle = 0 #Is the azimuthal angle left/right; change in script to modify for ears separately
alpha = 0 # in radians; useful if the movement is rotational
Vsound = 343 # Velocity of sound


Amp_file = './ref_table/Amp.txt'
yampl = np.loadtxt(Amp_file, delimiter='\t')

Freq_file = './ref_table/freqaxis.txt'
freqaxis = np.loadtxt(Freq_file, delimiter='\t')

Phase_file = './ref_table/yphase.txt'
yphase = np.loadtxt(Phase_file, delimiter='\t')
# print(yphase)
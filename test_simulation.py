from unittest import TestCase

from matplotlib import pyplot as plt


class Test(TestCase):
    def test_dist_to_surface(self):
        self.fail()

    def test_dist_to_surface_plot(self, source=None, sel_centroids=None, centroids=None, sel_centroids_2=None):
        # self.fail()
        print(sel_centroids.shape, centroids.shape)
        # checking the selected triangles plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection= '3d')

        ax.scatter(centroids[:,0], centroids[:,1], centroids[:,2], c='gray')
        ax.scatter(sel_centroids_2[:,0], sel_centroids_2[:,1], sel_centroids_2[:,2]+0.1, marker='*')
        ax.scatter(source[0], source[1], source[2], c = 'red')
        plt.show()



if __name__ == '__main__':
    from simulation import *
    # Configuration start
    Bat_to_Landscape = 8.5  # ditance to the landscape
    Azimut = 0  # head direction (depend on bats)
    elevation_angle = 90  # target angle to the center point = 180 - elevation_angle
    mikes_dist_coef = 0.0035  # pip kuli
    # Azimuth = 0

    # echo part
    min_freq = 35e3
    max_freq = 90e3  # for Pipistrellus kuhlii
    step_size = 50  # frequency step size
    freqvec = np.arange(min_freq, max_freq + 1, step_size)  # Emitter frequencies 25 kHz to 80 kHz
    # print(freqvec)

    Fs = 4096 * step_size * 2  # sample rate

    if Fs < 3 * max_freq:
        print('Warning, Raise step size')

    Em_level = 1  # use 1 Pa as emission
    # freq = 10 # has no function except for bemfa code
    max_gain_of_ear_db = 21
    mic_angle = 0  # Is the azimuthal angle left/right; change in script to modify for ears separately
    alpha = 0  # in radians; useful if the movement is rotational
    Vsound = 343  # Velocity of sound

    Amp_file = './ref_table/Amp.txt'
    yampl = np.loadtxt(Amp_file, delimiter='\t')

    Freq_file = './ref_table/freqaxis.txt'
    freqaxis = np.loadtxt(Freq_file, delimiter='\t')

    Phase_file = './ref_table/yphase.txt'
    yphase = np.loadtxt(Phase_file, delimiter='\t')
    # print(yphase)
    # Configuration end


    file = './test_obj/remeshed2mm.obj'
    # file = './test_obj/box.obj'
    # file = './test_obj/sphere.obj'

    # read the file
    vertices, faces_v = readObj(file)

    # center the vertices
    vertices = centerGeom(vertices)
    # print(vertices)

    # compute the normal using vertices and faces
    normals = getNormals(vertices, faces_v)

    centroids = getCentroids(vertices, faces_v)
    # print(centroids)

    # compute the average edge length
    triangle_edge_distances = SeparateEdgeLength(vertices, faces_v)
    # print(mean_edge_distance)

    # source information
    # print(x, y, z)
    # x, y, z = sph2cart(np.deg2rad(90 - Azimut), np.deg2rad(elevation_angle), Bat_to_Landscape)
    x, y, z = 0,0,8.5
    print(f'source: {x}, {y}, {z}')

    # distance of different ears of emiters because bats has two ears, which get somewhat different signals
    mike1, mike2 = calculate_ear_poisiton(x, y, z, Azimut, mikes_dist_coef)
    mic = calculate_reciver_poisiton(x, y, z, Azimut, mikes_dist_coef, number_of_mics=3)
    print(f'mic: {mic}')

    # to select triangles
    source = [x, y, z]
    sel_centroids, sel_projpoints, sel_source2centroids, fake_position = distance_to_surface(source, centroids, normals)
    print('selected triange number: ', len(sel_centroids), 'out of ', len(centroids))
    #
    # # calculate the distance:
    angles_mike1, angles_mike2, all_distancesTOWARDS, all_distancesbackmike1, all_distancesbackmike2 = \
        get_distances_and_angles(mike1, mike2, sel_centroids, sel_projpoints, sel_source2centroids)

    # calculate recivers distance:
    mics_angles, all_distancesTOWARDS, all_distances_back_mics = \
        get_reciver_distances_and_angles(mic, sel_centroids, sel_projpoints, sel_source2centroids)
    #
    # print("mics_angles: ", mics_angles)
    # print("all_distancesTOWARDS: ", all_distancesTOWARDS)
    # print('all_distancesbackmike1: ', all_distancesbackmike1)
    # # sys.exit(1)
    #
    # # select sonar beam range using HRTF
    HRTF_corr1, HRTF_corr2 = HRTFscript(alpha, mikes_dist_coef, Bat_to_Landscape, mic_angle, mike1, mike2, sel_centroids, freqvec, max_gain_of_ear_db)
    #
    # # echo_simulation
    # W3_left, W3_right, w3_left, w3_right = create_echo_single_triangles(triangle_edge_distances[fake_position],
    #                                                                     angles_mike1, angles_mike2,
    #                                                                     all_distancesTOWARDS, all_distancesbackmike1,
    #                                                                     all_distancesbackmike2, HRTF_corr1, HRTF_corr2)

pass
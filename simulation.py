# import sys
import numpy as np
from numpy import ndarray
import numpy.ma as ma
from matplotlib import pyplot as plt
import sys
from math import acos
from scipy.special import jn
import math

speed_of_sound = 343


# open the finite element in obj format file and organize the data
def yield_file(in_file):
    f = open(in_file)
    buf = f.read()
    f.close()

    # print(buf)
    for b in buf.split('\n'):
        # print('b: ', b)
        # sys.exit(1)
        if b.startswith('v '):
            yield ['v', [float(x) for x in b.split(" ")[1:4]]]
        elif b.startswith('f '):
            triangle = b.split(' ')[1:4]

            try:
                array_tmp = np.array([[int(j) - 1 for j in t.split("//")] for t in triangle if t != ''])
                array_tmp = array_tmp.T.reshape(2, 3)[0]
                # print('I am here!! 1111')
            except:
                array_tmp = np.array([int(t) - 1 for t in triangle if t != ''])
                # print('I am here!! 2222')

                # print('Choose alternative way to parse the obj file!!')
            # print(array_tmp)

            yield ['f', array_tmp]

        # elif b.startswith('vn '):
        # 	yield ['vn', [float(x) for x in b.split(" ")[1:]]]

        else:
            yield ['', ""]
    # NOTICE: don't return any object


# in Obj file,
# v: Vertex in format (x,y,z) (require)
# vn: Vertex normal in format (x,y,z) (optional)
# vt: Texture coordinates (UV Mapping) in format (u,v)
# projecting a 3D model's surface to a 2D image for texture mapping
# f: Face element (Polygonal) usually as triangle
# in format of "f x1 x2 x3 ..." when x can be v, v/vt, v1/vt1/vn1, v1//vn1
# it is written as numbers in the list e.g 1/2/5 v number 1,vt number 2, vn number 5
# take the name of the file and return array of vertical and face vector
# method: for loop
def readObj(objFileName):
    vertices = []
    faces_v = []
    faces_vt = []
    normals = []

    for k, v in yield_file(objFileName):
        if k == 'v':
            vertices.append(v)
        elif k == 'f':
            faces_v.append(list(v))
            # faces_vt.append(list(v[1]))
        # elif k == 'vn':
        # 	normals.append(v)

    # print(vertices)
    return np.array(vertices), np.array(faces_v)


# centering all the verticals vectors
# method: mean of the verices and
def centerGeom(vertices):
    # print(vertices)
    geomCentroid = np.mean(vertices, axis=0)
    cVertices = vertices - geomCentroid

    return cVertices


# take the vertices and faces and return normal vector array
# method:
def getNormals(vertices, faces):
    face_normals: ndarray = np.cross(vertices[faces[:, 1]] - vertices[faces[:, 0]],
                                     vertices[faces[:, 2]] - vertices[faces[:, 1]])
    face_areas = np.sqrt((face_normals ** 2).sum(axis=1))
    face_normals /= face_areas[:, np.newaxis]
    return face_normals


# centroid: the geometric center of the triangle
# in the triangle it is the Median
# NOTICE: calculation by the mean, it's not the write way
def getCentroids(vertices, faces):
    x = vertices[faces[:]].T.transpose(0, 2, 1)
    faceCentroids = np.mean(x, axis=2).T
    # print(x.shape, faceCentroids.shape)
    return faceCentroids


# take position and calculate the length
def edge_length(positions):
    return np.sqrt((positions ** 2).sum(axis=1))


# NOTICE: treat the edges like a triangle (assumption)
# take the
def SeparateEdgeLength(vertices, faces):
    tempVect1 = vertices[faces[:, 0]]
    tempVect2 = vertices[faces[:, 1]]
    tempVect3 = vertices[faces[:, 2]]

    edge1_length = edge_length(tempVect1 - tempVect2)
    edge2_length = edge_length(tempVect1 - tempVect3)
    edge3_length = edge_length(tempVect2 - tempVect3)
    tempVect_v = np.vstack((edge1_length, edge2_length, edge3_length)).T

    return np.mean(tempVect_v, axis=1)


# distance between centroids to to the sourse
def centroids_to_source_distance(source, centroids):
    # distances = []
    distances = [np.sqrt(np.sum((p - source) ** 2, axis=0)) for p in centroids]
    return np.array(distances)


# Spherical to Cartesian coordination
def sph2cart(azimuth, elevation, r):
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z


# take the source, azimut, and distance coeficient and number of reciver/2 and create a recivers array locations
def calculate_reciver_poisiton(x, y, z, Azimut, distance_to_mics, number_of_mics=1):
    # n - will be duplicated to the two side
    array_center_location = [x, y, z]
    mic = []
    for i in range(number_of_mics):
        sign = +1 if ((i + 1) % 2) == 0 else -1
        mic.append([x + sign * np.cos(np.deg2rad(Azimut)) * distance_to_mics * (i // 2 + 1),
                    y - sign * np.sin(np.deg2rad(Azimut)) * distance_to_mics * (i // 2 + 1), z])
    return mic

def calculate_reciver_poisiton_2d(x, y, z, Azimut, distance_to_mics, number_of_mics=1):
    # n - will be duplicated to the two side
    array_center_location = [x, y, z]
    mic = []
    for i in range(number_of_mics):
        sign = +1 if ((i + 1) % 2) == 0 else -1
        mic.append([x + sign * np.cos(np.deg2rad(Azimut)) * distance_to_mics * (i // 2 + 1),
                    y - sign * np.sin(np.deg2rad(Azimut)) * distance_to_mics * (i // 2 + 1), z])
    return mic

def level_def(segmlen):
    global level
    if segmlen > 0.002 and segmlen < 0.02:
        level = 12947 * segmlen + 51.967

    if segmlen > 0.02:
        level = 280

    if segmlen < 0.002:
        level = 0.0081 * segmlen ** (-1.482)

    return level


# finding all the triangle that projected from the source
# Notice: to be more accurate needs to include all the sources / feadback from other object
def distance_to_surface(source, centroids, normals):
    # start = time.time()
    source2centroids = centroids_to_source_distance(source, centroids)
    # print(source2centroids.shape, normals.shape)

    proj_points = centroids - (source2centroids * normals.T).T

    # checking the position
    z_mean_position = np.nanmean(proj_points[:, 2])
    if z_mean_position < 0:
        # print(z_mean_position)
        print('Notice: the projrction points is in minus, it is fixed .')
        proj_points = centroids + (source2centroids * normals.T).T

    delta_vecsq = (abs(centroids - proj_points)) ** 2
    centroids2projpoints = np.sqrt(sum(delta_vecsq.T))
    source2projpoints = centroids_to_source_distance(source, proj_points)

    hypot = np.sqrt(source2centroids ** 2 + centroids2projpoints ** 2)
    # ideal     # actual
    selected_indexes = np.argwhere(hypot > 1.05 * source2projpoints).T[0]

    # sel_normals = normals[selected_indexes]
    sel_centroids = centroids[selected_indexes]
    sel_projpoints = proj_points[selected_indexes]
    sel_source2centroids = source2centroids[selected_indexes]

    # print('selected triangle number: ', len(selected_indexes))
    if len(selected_indexes) < len(centroids) / 6:
        print('checking the selected the triangles because the selected triangles are too low', len(selected_indexes))
        sys.exit(1)

    # how to find occluder-object triangles??
    longest_objectaxis = max(source2centroids) - min(source2centroids)

    n = 4  # for scenes with a long depth you need to increase this number (4 suited for 1 object)

    potential_value = len(source2centroids)
    segmlen = (longest_objectaxis) / n
    level = level_def(segmlen)
    temp_source2centroids = source2centroids

    centr_ind = np.array([np.nan] * potential_value * n).reshape(n, potential_value)
    # print(centr_ind)

    # sorted_s2c = sorted(source2centroids)
    my_segment_numbers = [sum(source2centroids > (max(source2centroids) - segmlen * (i + 1) * 1.001)) for i in range(n)]
    sorted_s2c_index = np.flipud(np.argsort(source2centroids))  # get index and order from max to min

    selraddisected_objectsegm = []
    allraddisected_objectsegm = []
    for idx, number in enumerate(my_segment_numbers):
        if idx == 0:
            tmp = list(sorted(sorted_s2c_index[:number]))
        else:
            tmp = list(sorted(sorted_s2c_index[my_segment_numbers[idx - 1]:my_segment_numbers[idx]]))

        allraddisected_objectsegm.append(centroids[tmp])

        filter_preind = sorted(list(set(selected_indexes) & set(tmp)))
        selraddisected_objectsegm.append(centroids[filter_preind])

        # print('1111: ', filter_preind)
        # sys.exit(1)

        for i in range(len(filter_preind)):
            centr_ind[idx][i] = filter_preind[i]

    # this is core of this algorithm - to select the potential triangle without occludes
    Occl_vec = []
    for pik in range(n - 1):
        for kut in range(pik + 1, n):
            emitting_faces = selraddisected_objectsegm[pik]
            occluders = allraddisected_objectsegm[kut]

            direc = source - emitting_faces  # direction to the emitting faces

            ocl_x = occluders[:, 0]
            ocl_y = occluders[:, 1]
            ocl_z = occluders[:, 2]

            # x = (np.array([emitting_faces[:,0]]*len(ocl_x.T)).T - np.array([ocl_x.T]*len(emitting_faces[:,0])))*10000
            t_x = (ocl_x.T - emitting_faces[:, 0][:, None]) * 10000 / direc[:, 0][:, None]
            t_y = (ocl_y.T - emitting_faces[:, 1][:, None]) * 10000 / direc[:, 1][:, None]
            t_z = (ocl_z.T - emitting_faces[:, 2][:, None]) * 10000 / direc[:, 2][:, None]

            t_x = ma.masked_where((t_x > -50) & (t_x < 50), t_x)
            t_y = ma.masked_where((t_y > -50) & (t_y < 50), t_y)
            t_z = ma.masked_where((t_z > -50) & (t_z < 50), t_z)
            posxy = ma.masked_greater_equal(abs(t_x - t_y), level).mask + np.ma.masked_greater_equal(abs(t_x - t_z),
                                                                                                     level).mask
            idxs = np.argwhere(np.any(~posxy, axis=1) == True)

            if len(idxs) > 0:
                occlIndeces = centr_ind[pik, idxs].ravel()
                Occl_vec.extend([int(i) for i in occlIndeces])

    # print('Occl_vec: ', len(Occl_vec))
    fake_position = sorted(list(set(selected_indexes) - set(Occl_vec)))
    # print(fake_position, len(fake_position))

    # sel_centroids_2 = centroids[fake_position]
    # sel_normals_2 = normals[fake_position]
    sel_centroids_2 = centroids[fake_position]
    sel_projpoints_2 = proj_points[fake_position]
    sel_source2centroids_2 = source2centroids[fake_position]
    return sel_centroids_2, sel_projpoints_2, sel_source2centroids_2, fake_position


def check_angle(angles):
    my_min = min(angles)
    my_max = max(angles)
    if my_min < 0 or my_max > 91:
        angles[angles > 90] = 90
        print('min angle: ', my_min, '\tmax angle:', my_max)
        print('Warning!!!The angle is wrong, please check the selection of the triangle!!')
    # sys.exit(1)


def get_reciver_distances_and_angles(mic_array, sel_centroids, sel_projpoints, sel_source2centroids):
    global all_distances_forward
    delta_vecsq = (abs(sel_centroids - sel_projpoints)) ** 2
    a_c2p = np.sqrt(sum(delta_vecsq.T))
    mics_angles = []
    all_distances_back_mics = []
    for i, mic in enumerate(mic_array):
        b_c2mic = centroids_to_source_distance(mic, sel_centroids)
        c_p2mic = centroids_to_source_distance(mic, sel_projpoints)

        tmp = (a_c2p ** 2 + b_c2mic ** 2 - c_p2mic ** 2) / (2 * a_c2p * b_c2mic)
        mics_angles.append(np.rad2deg([acos(j) for j in tmp]))

        all_distances_forward = sel_source2centroids
        all_distances_back_mics.append(b_c2mic.T)

        check_angle(mics_angles[-1])

        # angles_mike1[angles_mike1 > 1] = 90

    return mics_angles, all_distances_forward, all_distances_back_mics


def HRTF(Foc_Point_mike, mic, min_scale_highest_frec, a, waven, sel_centroids, max_gain_of_ear_db):
    # Calculate distances between microphone and focal point, microphone and centroids, and line connecting them
    ov = centroids_to_source_distance(Foc_Point_mike, sel_centroids)
    schuine = centroids_to_source_distance(mic, sel_centroids)
    linec = centroids_to_source_distance(mic, Foc_Point_mike)[1]

    # Compute cosine of the angle between microphone-centroid line and focal point-centroid line
    cosine_value = ((schuine ** 2 - ov ** 2) + linec ** 2) / (2 * schuine * linec)

    # Avoid invalid values for arccosine (due to floating-point precision)
    cosine_value[cosine_value > 1] = 1

    # Compute angle in degrees relative to main-axis to Foc-point
    angle_rad = np.arccos(cosine_value)
    angle_deg = np.degrees(angle_rad)
    angle_deg[angle_deg > 180] = 180  # Restrict angles to faces in front

    # Calculate HRTF corrections
    tmp = waven[:, np.newaxis] * np.sin(np.radians(angle_deg) + 1e-9) * a
    HRTF_corr = max_gain_of_ear_db * (
            (2 * jn(1, tmp) / tmp + abs(min_scale_highest_frec)) / (1 + abs(min_scale_highest_frec)))

    return HRTF_corr


def get_position(freqvec, freqaxis, edge_distance):
    return np.array(
        [np.argmin(abs(freq_tmp[:, None] - freqaxis), axis=1) for freq_tmp in freqvec[:, None] * edge_distance / 0.1])


# correction of the recivers (Need to be chacked im not sure about this function)
def HRTF_receiver_script(alpha, mic_distance_coefficient, ideal_dist_to_obj, mic_angle, mic_array, sel_centroids,
                         freq_vec, max_gain_of_ear_db, emitter_radius=0.006):
    speed_of_sound = 343
    wavenumber = 2 * np.pi / (speed_of_sound / freq_vec)

    tmp = wavenumber[-1] * emitter_radius * np.sin(np.deg2rad(np.arange(-90, 91)))  # 0 to 180 degrees
    scaling_bessel = 2 * jn(1, tmp) / tmp

    min_scale_highest_freq = min(scaling_bessel)

    Foc_Point_source = [0, 0, 0]
    HRTF_corr_array = []

    for i, mic_position in enumerate(mic_array):
        sign = +1 if (i + 1) % 2 == 0 else -1
        foc_distance = np.tan(np.deg2rad(mic_angle)) * ideal_dist_to_obj
        Foc_Point_mic = [
            Foc_Point_source[0] + sign * np.cos(alpha) * (mic_distance_coefficient * (i // 2 + 1) + foc_distance),
            Foc_Point_source[1] + sign * np.sin(alpha) * (mic_distance_coefficient * (i // 2 + 1) + foc_distance),
            0
        ]
        HRTF_corr_array.append(
            HRTF(np.array(Foc_Point_mic), np.array(mic_position), min_scale_highest_freq, emitter_radius, wavenumber,
                 sel_centroids, max_gain_of_ear_db))

    return np.array(HRTF_corr_array)


def HRTF_omni(mic_position, min_scale_highest_frequncy, emitter_radius, sel_centroids, max_gain_of_ear_db):
    # For an omnidirectional microphone, no HRTF correction needed
    # You can simply return the maximum gain for all frequencies

    # Calculate distances from mic to centroids
    distances = centroids_to_source_distance(mic_position, sel_centroids)

    # Calculate HRTF correction
    tmp = distances[1] * emitter_radius
    HRTF_corr = max_gain_of_ear_db * (
            (2 * jn(1, tmp) / tmp + abs(min_scale_highest_frequncy)) / (1 + abs(min_scale_highest_frequncy)))
    return HRTF_corr


def HRTF_receiver_script_omni(mic_array, sel_centroids, freq_vec, max_gain_of_ear_db, emitter_radius=0.006):
    speed_of_sound = 343
    # Calculate wavenumber for each frequency
    wavenumber = 2 * np.pi * freq_vec / speed_of_sound

    # Calculate Bessel function for scaling
    tmp = wavenumber[-1] * emitter_radius * np.arange(0,
                                                      len(sel_centroids))  # Assuming omnidirectional, angles don't matter

    # Handle division by zero or small values in tmp
    mask = np.isclose(tmp, 0)  # Find elements close to zero
    scaling_bessel = np.zeros_like(tmp)  # Initialize scaling_bessel with zeros
    scaling_bessel[~mask] = 2 * jn(1, tmp[~mask]) / tmp[
        ~mask]  # Calculate scaling_bessel where tmp is not close to zero
    # For elements close to zero, scaling_bessel will remain zero

    # Find minimum scaling factor at the highest frequency
    min_scale_highest_freq = np.nanmin(scaling_bessel)  # Ignore NaN values while finding minimum

    HRTF_corr_array = []

    # Calculate HRTF correction for each omnidirectional microphone
    for mic_position in mic_array:
        HRTF_corr_array.append(HRTF_omni(np.array(mic_position), min_scale_highest_freq, emitter_radius, sel_centroids,
                                         max_gain_of_ear_db))

    return np.array(HRTF_corr_array)


def get_values_from_table(table, positions, angles):
    return np.array([table[(angles, position)] for position in positions])  # 321*9518


# def create_echo_single_triangles(edgeLengths, angles_mikes, all_distancesTOWARDS, all_distancesbackmikes,HRTF_corr):
#     # calculate total reflection spectrum due to all faces
#     W = Vsound / freqvec
#
#     # Edge lengths and distances
#     edgeLengths = edgeLengths
#     takethisposs = get_position(freqvec, freqaxis, edgeLengths)
#
#     combloss_dB = []
#     ang_shifts = []
#     losses = []
#
#     # Calculate loss and phase shift for each receiver
#     for angles_mike in angles_mikes:
#         angles_mike = np.array([round(abs(i)) for i in angles_mike])
#         combloss_dB.append(20 * np.log10(get_values2(yampl, takethisposs, angles_mike)))
#         ang_shifts.append(get_values2(yphase, takethisposs, angles_mike))
#
#     # Calculate total distances for each receiver
#     tot_dists = all_distancesTOWARDS + all_distancesbackmikes
#     tot_dists = ma.masked_where(tot_dists < 0, tot_dists)
#
#     # Calculate phase and atmospheric losses for each receiver
#     phase_shifts = []
#     atm_losses = []
#     for tot_dist in tot_dists:
#         tot_dist_mat, wavele_mat = np.meshgrid(tot_dist, W)
#         fract_w = tot_dist_mat / wavele_mat
#         phase = (fract_w % 1) * 2 * math.pi
#         phase_shifts.append(phase)
#         atm_losses.append(-2.591393338e-25 * freqvec ** 5 + 1.040028037e-19 * freqvec ** 4 -
#                           1.46257453e-14 * freqvec ** 3 + 8.652366167e-10 * freqvec ** 2 +
#                           1.627757516e-5 * freqvec - 0.11795019 * tot_dist)
#
#     # Calculate total losses for each receiver
#     for i in range(len(angles_mikes)):
#         loss = 2 * 6.02 * np.log2(0.08 / (tot_dists[i] / 2))
#         losses.append(loss + combloss_dB[i] + HRTF_corr - atm_losses[i])
#
#     # Calculate phase of each arriving ray for each receiver
#     returns = []
#     for i in range(len(angles_mikes)):
#         phase_b = np.exp(-1j * (math.pi + phase_shifts[i] + ang_shifts[i]))
#         array = np.zeros((len(tot_dists[i]), 1)) + Em_level
#         returns.append(phase_b * array.T)
#
#     # Convert losses to Pascal
#     outvolts = []
#     for i in range(len(angles_mikes)):
#         outvolt = returns[i] * (10 ** (losses[i] / 20))
#         outvolts.append(outvolt)
#
#     # Sum pressures for each receiver
#     strengths = [np.nansum(outvolt, 1) for outvolt in outvolts]
#
#     # Create spectrum impulse response for each receiver
#     n2 = 4096
#     Powers = [np.zeros(n2, dtype="complex_") for _ in range(len(angles_mikes))]
#     for i in range(len(angles_mikes)):
#         if len(freqvec) == len(strengths[i]):
#             indexes = [int(j - 1) for j in freqvec / step_size]
#             Powers[i][indexes] = strengths[i]
#         else:
#             print('the length of freqvec and strength is not equal')
#             sys.exit(1)
#
#     Hs = []
#     for Power in Powers:
#         H = np.flipud(np.conj(Power[1:n2]))
#         Hs.append(H)
#
#     W3s = []
#     w3s = []
#     for H in Hs:
#         W3 = np.hstack((Power, [0 + 0j], H))
#         w3 = np.real(np.array([np.fft.ifft(W3, n=len(W3))]).T)
#         W3s.append(W3)
#         w3s.append(w3)
#
#     return W3s, w3s


def create_echo_single_triangles(edgeLengths, mics_angles_array, all_distancesTOWARDS, all_distances_back_mics,
                                 HRTF_corr_array, source_frequancy_vector, parameters):
    # Calculate total reflection spectrum due to all faces
    W = speed_of_sound / source_frequancy_vector
    wave_power = []
    sound_wave = []
    yampl, freqaxis, yphase, Em_level, fft_size, step_size = parameters
    for i, mic_angles in enumerate(mics_angles_array):
        mic_angles = np.array([round(abs(i)) for i in mic_angles])
        takethisposs = get_position(source_frequancy_vector, freqaxis, edgeLengths)
        comblossdB = 20 * np.log10(get_values_from_table(yampl, takethisposs, mic_angles))
        ang_shift = get_values_from_table(yphase, takethisposs, mic_angles)

        tot_dist = all_distancesTOWARDS + all_distances_back_mics[i]
        tot_dist = ma.masked_where(tot_dist < 0, tot_dist)
        tot_dist_mat, wavele_mat = np.meshgrid(tot_dist, W)

        # Create frequency-dependent atmospheric losses
        Atm = -2.591393338e-25 * source_frequancy_vector ** 5 + 1.040028037e-19 * source_frequancy_vector ** 4 - 1.46257453e-14 * source_frequancy_vector ** 3 + \
              8.652366167e-10 * source_frequancy_vector ** 2 + 1.627757516e-5 * source_frequancy_vector - 0.11795019
        Atm_loss = Atm[:, np.newaxis] * tot_dist
        loss = 2 * 6.02 * np.log2(0.08 / (tot_dist / 2))
        Ldir_path = (loss + comblossdB + HRTF_corr_array[i]) - Atm_loss

        # Calculate phase of each arriving ray for each receiver
        fract_w = tot_dist_mat / wavele_mat
        phase = (fract_w % 1) * 2 * math.pi
        phaseb = np.exp(-1j * (math.pi + phase + ang_shift))
        array = np.zeros((len(tot_dist), 1)) + Em_level
        returns = phaseb * array.T

        # Convert losses to Pascal
        outvolt = returns * (10 ** (Ldir_path / 20))
        # Sum pressures for each receiver
        strength_ear = np.nansum(outvolt, 1)
        # Create spectrum impulse response for each receiver
        # fft_size = 4096
        Power = np.zeros(fft_size, dtype="complex_")
        if len(source_frequancy_vector) == len(strength_ear):
            indexes = [int(i - 1) for i in source_frequancy_vector / step_size]
            Power[indexes] = strength_ear
        else:
            print('the length of freqvec and strength is not equal')
            sys.exit(1)

        # this action is not understood
        H = np.flipud(np.conj(Power[1:fft_size]))

        wave_power.append(np.hstack((Power, [0 + 0j], H)))
        sound_wave.append(np.real(np.array([np.fft.ifft(np.array(wave_power[-1]), n=len(wave_power[-1]))]).T))

    return np.array(wave_power), sound_wave


def create_echo_s_triangles(edgeLengths, mic_angles, all_distancesTOWARDS, all_distances_back_mics, HRTF_corr_array,
                            source_frequancy_vector,
                            source_frequency_vector=None):
    # Calculate total reflection spectrum due to all faces
    W = speed_of_sound / source_frequancy_vector
    wave_power = []
    sound_wave = []
    for i, mic_angles in enumerate(mics_angles_array):
        mic_angles = np.array([round(abs(i)) for i in mic_angles])
        takethisposs = get_position(source_frequancy_vector, freqaxis, edgeLengths)
        comblossdB = 20 * np.log10(get_values_from_table(yampl, takethisposs, mic_angles))
        ang_shift = get_values_from_table(yphase, takethisposs, mic_angles)

        tot_dist = all_distancesTOWARDS + all_distances_back_mics[i]
        tot_dist = ma.masked_where(tot_dist < 0, tot_dist)
        tot_dist_mat, wavele_mat = np.meshgrid(tot_dist, W)

        # Create frequency-dependent atmospheric losses
        Atm = -2.591393338e-25 * source_frequancy_vector ** 5 + 1.040028037e-19 * source_frequancy_vector ** 4 - 1.46257453e-14 * source_frequancy_vector ** 3 + \
              8.652366167e-10 * source_frequancy_vector ** 2 + 1.627757516e-5 * source_frequancy_vector - 0.11795019
        Atm_loss = Atm[:, np.newaxis] * tot_dist
        loss = 2 * 6.02 * np.log2(0.08 / (tot_dist / 2))
        Ldir_path = (loss + comblossdB + HRTF_corr_array[i]) - Atm_loss

        # Calculate phase of each arriving ray for each receiver
        fract_w = tot_dist_mat / wavele_mat
        phase = (fract_w % 1) * 2 * math.pi
        phaseb = np.exp(-1j * (math.pi + phase + ang_shift))
        array = np.zeros((len(tot_dist), 1)) + Em_level
        returns = phaseb * array.T

        # Convert losses to Pascal
        outvolt = returns * (10 ** (Ldir_path / 20))
        # Sum pressures for each receiver
        strength_ear = np.nansum(outvolt, 1)
        # Create spectrum impulse response for each receiver
        # fft_size = 4096
        Power = np.zeros(fft_size, dtype="complex_")
        if len(source_frequency_vector) == len(strength_ear):
            indexes = [(int(freq / step_size) - 1) for freq in source_frequency_vector]
            Power[indexes] = strength_ear
        else:
            print('The length of source_frequency_vector and strength is not equal')
            sys.exit(1)
        # this action is not understood
        H = np.flipud(np.conj(Power[1:fft_size]))

        wave_power.append(np.hstack((Power, [0 + 0j], H)))
        sound_wave.append(np.real(np.array([np.fft.ifft(np.array(wave_power[-1]), n=len(wave_power[-1]))]).T))

    return np.array(wave_power), sound_wave


def create_echoes(waves, sample_rate, echo_intervals):
    """
    Create echoes of sound waves for multiple receivers with specified intervals.
    """
    combined_waves_with_echoes = []

    for wave in waves:
        # Create empty list to store delayed echoes for this receiver
        echoes = []
        max_echo_time = int(np.max(echo_intervals) * sample_rate + len(wave))
        # Apply delays to the original sound wave for this receiver
        for interval in echo_intervals:
            # Calculate the number of samples to delay based on the delay interval and sample rate
            delay_samples = int(interval * sample_rate)
            # Create a delayed version of the original sound wave using zero-padding
            # delayed_wave = delayed_empty
            delayed_wave = np.concatenate((np.zeros(delay_samples), np.array(wave.flatten()),
                                           np.zeros(max_echo_time - delay_samples - len(wave))))

            # Append the delayed version to the list of echoes
            echoes.append(delayed_wave)

        # Combine the original sound wave with the echoes for this receiver
        combined_wave = np.vstack(np.array(echoes)).sum(axis=0)

        # Append the combined wave with echoes to the list
        combined_waves_with_echoes.append(combined_wave)

    return combined_waves_with_echoes


if __name__ == '__main__':

    # Configuration start

    # OBJECT
    # file = './test_obj/remeshed2mm.obj'
    # file = './test_obj/sample.obj'
    # file = './test_obj/box.obj'
    file = './test_obj/sphere.obj'

    # Recivers Setting
    reciver_to_object = 8.5  # ditance to the landscape
    azimut = 0  # head direction (depend on bats)
    elevation_angle = 90  # target angle to the center point = 180 - elevation_angle
    distance_between_mics = 1  # pip kuli
    # distance_between_mics = 0.035  # pip kuli
    number_of_mics = 8

    # source location
    x_source, y_source, z_source = 0, 0, reciver_to_object
    source = [x_source, y_source, z_source]
    print(f'source: {x_source}, {y_source}, {z_source}')

    # recive location
    # distance of different ears of emiters because bats has two ears, which get somewhat different signals
    mic_array = calculate_reciver_poisiton(x_source, y_source, z_source, azimut, distance_between_mics, number_of_mics)
    print(f'mic: {mic_array}')

    # Source Setting
    min_frequncy = 25e3
    max_frequncy = 60e3  # for Pipistrellus kuhlii
    step_size = 50  # frequency step size
    source_frequancy_vector = np.arange(min_frequncy, max_frequncy + 1,
                                        step_size)  # Emitter frequencies 25 kHz to 80 kHz

    # simulation
    fft_step_size = 64  # frequency step size 2*6
    fft_size = 4096
    sample_rate = fft_size * fft_step_size  # sample rate (multiply in 2)
    # sample_rate = 44100  # sample rate
    time_interval = 1 / sample_rate  # sec
    # total_time = fft_size * time_interval
    total_time = 1 / fft_step_size
    print(f"total time simulation: {total_time}")
    window_size = 512  # Window size for FFT (2^9)
    overlap = 256  # Overlap for FFT (2^8)

    mic_angle = 0  # Is the azimuthal angle left/right; change in script to modify for ears separately
    alpha = 0  # in radians; useful if the movement is rotational
    max_gain_of_ear_db = 21
    Em_level = 1
    speed_of_sound = 343

    # Warning of the setting
    if sample_rate < 3 * max_frequncy:
        print('Warning, Raise step size')

    # Sound pressure

    # print(source_frequancy_vector)

    # referance tables
    Amp_file = './ref_table/Amp.txt'
    yampl = np.loadtxt(Amp_file, delimiter='\t')

    Freq_file = './ref_table/freqaxis.txt'
    freqaxis = np.loadtxt(Freq_file, delimiter='\t')

    Phase_file = './ref_table/yphase.txt'
    yphase = np.loadtxt(Phase_file, delimiter='\t')

    parameters = [yampl,freqaxis,yphase,Em_level,fft_size,step_size]
    # print(yphase)
    # Configuration end

    # read the file
    vertices, faces_v = readObj(file)

    # center the vertices
    vertices = centerGeom(vertices)
    # print(vertices)

    # compute the normal using vertices and faces
    normals = getNormals(vertices, faces_v)  # normals of the object Finate Element tringle
    centroids = getCentroids(vertices, faces_v)  # centroid of each tringle
    # print(centroids)

    # compute the average edge length
    triangle_edge_distances = SeparateEdgeLength(vertices, faces_v)
    # print(mean_edge_distance)

    # select all the triangle that projected from the source
    sel_centroids, sel_projpoints, sel_source2centroids, fake_position = \
        distance_to_surface(source, centroids, normals)
    print('selected triange number: ', len(sel_centroids), 'out of ', len(centroids))

    # calculate recivers distance:
    mics_angles_array, all_distancesTOWARDS_r, all_distances_back_mics = get_reciver_distances_and_angles(mic_array,
                                                                                                          sel_centroids,
                                                                                                          sel_projpoints,
                                                                                                          sel_source2centroids)

    HRTF_omni_array = HRTF_receiver_script_omni(mic_array, sel_centroids, source_frequancy_vector, max_gain_of_ear_db)

    #
    # # Sonar_simulation with single
    wave_power_array, sound_wave_array = create_echo_single_triangles(triangle_edge_distances[fake_position],
                                                                      mics_angles_array, all_distancesTOWARDS_r,
                                                                      all_distances_back_mics, HRTF_omni_array,
                                                                      source_frequancy_vector,parameters)
    # print(wave_power, wave_power.shape)

    # # plot of the sound wave from both microphone
    fig, ax = plt.subplots(len(sound_wave_array))
    fig.suptitle('Sound wave format')
    for i, sound_wave in enumerate(sound_wave_array):
        Time_axis = np.linspace(0, len(sound_wave) / 2 * sample_rate, len(sound_wave))
        ax[i].plot(Time_axis, sound_wave)
    # plt.show()
    #
    # Sonar_simulation with interval from each echo
    echo_interval = [0, 0.01,0.02,0.03,0.04]
    # echo_interval = [0, 0.05, 0.1, 0.15]
    combined_waves_with_echoes = create_echoes(sound_wave_array, sample_rate, echo_interval)
    # # plot of the sound wave from both microphone
    fig, ax = plt.subplots(len(combined_waves_with_echoes))
    fig.suptitle('Sound wave format')
    for i, ww in enumerate(combined_waves_with_echoes):
        Time_axis = np.linspace(0, len(ww) / sample_rate, len(ww))
        ax[i].plot(Time_axis, ww)
    # plt.show()

    # sample_rate = fft_size * fft_step_size  # Sample rate (multiply in 2)
    # window_size = 512  # Window size for FFT (2^9)
    # overlap = 256  # Overlap for FFT (2^8)
    #
    # # Plot sonogram for each receiver
    # for i, sound_wave in enumerate(sound_wave_array):
    #     plt.figure(figsize=(10, 6))  # Adjust figsize as needed for high resolution
    #     flattened_sound_wave = sound_wave.ravel()
    #     plt.specgram(flattened_sound_wave, Fs=sample_rate, NFFT=window_size, noverlap=overlap)
    #     plt.title(f'Sonogram - Receiver {i + 1}')
    #     plt.xlabel('Time (s)')
    #     plt.ylabel('Frequency (Hz)')
    #     plt.colorbar(label='Intensity (dB)')
    #     plt.tight_layout()
    #
    #     # Save the sonogram as a high-resolution image
    #     plt.savefig(f'./graph/receiver_{i + 1}_sonogram.png', dpi=300)  # Adjust dpi for high resolution

    plt.show()
    pass

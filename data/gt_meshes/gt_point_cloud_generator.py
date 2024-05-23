import trimesh
import numpy as np


def export_gt_point_cloud(mesh_file, export_name, gt_n_points = 20000, gt_n_points_volume = 50000):
    # Load G.T mesh
    gt_mesh = trimesh.load(mesh_file+".obj", force='mesh')
    gt_points, gt_faces = gt_mesh.sample(gt_n_points, return_index=True)
    gt_normals = gt_mesh.face_normals[gt_faces]
    np.save(export_name+"_surface", gt_points)
    np.save(export_name+"_surface_normal", gt_normals)

    gt_points_volume = trimesh.sample.volume_mesh(gt_mesh, gt_n_points_volume)
    np.save(export_name+"_volume", gt_points_volume)


#export_gt_point_cloud("bunny", 20000)
#export_gt_point_cloud("armadilo", 20000)
# export_gt_point_cloud("budda_very_low", "budda", 20000)
# export_gt_point_cloud("dragon_very_low", "dragon", 20000)
# export_gt_point_cloud("lucy_very_low", "lucy", 20000)
# export_gt_point_cloud("xyz_dragon_very_low", "xyz_dragon", 20000)
export_gt_point_cloud("cylinder", "cylinder", 20000)
export_gt_point_cloud("cube", "cube", 20000)

# export_gt_point_cloud("budda_low_res", "budda", 20000)
# export_gt_point_cloud("dragon_low_res", "dragon", 20000)
# export_gt_point_cloud("lucy_low_res", "lucy", 20000)
# export_gt_point_cloud("xyz_dragon_low_res", "xyz_dragon", 20000)

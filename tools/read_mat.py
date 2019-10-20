import scipy.io as sio
import numpy as np
import json
import os


def load_data(data_name):
    data = sio.loadmat(data_name)
    prediction = data['prediction'][0,0]
    image_path, image_size, line_segs, hvps, hgroup, zvps, zgroup, focal  = prediction


    image_path = image_path[0]
    image_size = image_size.tolist()
    image_size = [image_size[0][0], image_size[0][1]]  # height x width
    line_segs = line_segs.tolist()
    vps = np.concatenate((hvps, zvps), axis=0).tolist()
   
    group = -np.ones(len(line_segs)).astype(np.int)
    g_ind = 0
    for item in hgroup:
        h_list = item[0].T[0] - 1
        for h in h_list:
            group[h] = g_ind
        g_ind += 1
    zg = zgroup[0].T - 1
    for z in zg:
        group[z] = g_ind
    group = group.tolist()

    focal = float(focal[0][0])
    
    return image_path, image_size, line_segs, vps, group, focal


def point2line(end_points):
    # line: ax + by + c = 0, in which a^2 + b^2=1, c>0
    # point: 2 x 2  # point x dim
    # A = np.matrix(end_points) - np.array(image_size) / 2
    # result = np.linalg.inv(A) * np.matrix([1,1]).transpose()

    A = np.asmatrix(end_points)
    result = np.linalg.inv(A) * np.asmatrix([-1, -1]).transpose()  # a, b, 1
    a = float(result[0])
    b = float(result[1])
    norm = (a ** 2 + b ** 2) ** 0.5
    result = np.array([a / norm, b / norm, 1 / norm])

    return result


def lineseg2line(line_segs, image_size):
    # line_segs: number x (width, heigth)
    height, width = image_size
    new_line_segs = []
    new_lines = []
    for line_s in line_segs:
        end_points = [[line_s[1], line_s[0]], 
                      [line_s[3], line_s[2]]]
        new_line_segs.append(end_points)
        new_end_points = [[(end_points[i][0] - image_size[0] / 2 ) / (image_size[0] / 2),
                            (end_points[i][1] - image_size[1] / 2 ) / (image_size[1] / 2)]
                            for i in range(2)]
        new_line = point2line(new_end_points).tolist()
        new_lines.append(new_line)

    return new_line_segs, new_lines
        

def process(data_list, save_path):
    save_op = open(save_path, 'w')

    for data_name in data_list:
        print(data_name)
        image_path, image_size, line_segs, vps, group, focal = load_data(data_name)
       
        # image_size: height x width
        vps_output = []
        scale = max(image_size) / 2
        for vp in vps:
            new_vp = [(vp[1] - image_size[0] / 2) / (image_size[0] / 2), (vp[0] - image_size[1] / 2) / (image_size[1] / 2)]
            vps_output.append(new_vp)
        
        line_segs_output, new_lines_output = lineseg2line(line_segs, image_size)
        group_output = group
        
        image_name = image_path.split('/')[-1]
        json_out = {'image_path': image_name, 'line': new_lines_output, 'org_line': line_segs_output, 
                'group': group_output, 'vp': vps_output, 'focal': focal} 

        json.dump(json_out, save_op)
        save_op.write('\n')


if __name__ == '__main__':
    path = '/n/fs/vl/xg5/workspace/baseline/Horizon-First-VPdetection/dataset/YUD/output'
    data_list = [os.path.join(path, dir_path + '/data.mat') for dir_path in os.listdir(path)]

    save_path = '/n/fs/vl/xg5/workspace/baseline/Horizon-First-VPdetection/dataset/YUD/data/data.json'
    process(data_list, save_path)
    


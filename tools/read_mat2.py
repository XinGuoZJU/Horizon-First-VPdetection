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
    vps = np.concatenate((hvps, zvps), axis=0).tolist()
   
    focal = float(focal[0][0])
    
    return image_path, image_size, vps, focal


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
        image_path, image_size, vps, focal = load_data(data_name)
       
        # image_size: height x width
        vps_output = []
        scale = max(image_size) / 2
        for vp in vps:
            new_vp = [(vp[1] - image_size[0] / 2) / (image_size[0] / 2), (vp[0] - image_size[1] / 2) / (image_size[1] / 2)]
            vps_output.append(new_vp)
        
        image_names = image_path.split('/')
        image_name = os.path.join(image_names[-2], image_names[-1])
        
        json_out = {'image_path': image_name, 'image_size': image_size, 'vp': vps_output, 'focal': focal} 

        json.dump(json_out, save_op)
        save_op.write('\n')


if __name__ == '__main__':
    data_name = 'YUD'   # 'YUD', 'ScanNet', 'SceneCityUrban3D', 'SUNCG'

    path = '/n/fs/vl/xg5/workspace/baseline/Horizon_First_VPdetection/dataset/' + data_name + '/output'
    dir_list = [os.path.join(path, dir_path) for dir_path in os.listdir(path)]
    data_list = []
    for dirs in dir_list:
        data_list += [os.path.join(dirs, dir_path + '/data.mat') for dir_path in os.listdir(dirs)]

    save_path = '/n/fs/vl/xg5/workspace/baseline/Horizon_First_VPdetection/dataset/' + data_name + '/data'
    os.makedirs(save_path, exist_ok=True)
    save_file = os.path.join(save_path, 'data.json')

    process(data_list, save_file)
    


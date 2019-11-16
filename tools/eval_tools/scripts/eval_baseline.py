import os
import json
import numpy as np
from vp2focal2 import vp2focal


def nearest(p, points):
    # p: 2 (1D), points: number x 2
    distance = ((np.array(p) - np.array(points)) ** 2).sum(1).min() ** 0.5

    return float(distance)


def run(file_name, save_name, labels):
    label_dict = {}
    for label in labels:
        context = label.split()
        label_dict[context[0]] = [float(context[1]), [[float(context[2*i]), float(context[2*i+1])] for i in range(1,4)]]

    with open(file_name, 'r') as r_op:
        lines = r_op.readlines()
        vps_recall_list = []
        vps_precision_list = []
        focal_diff_list = []
        for line in lines:
            data_dict = json.loads(line)
            image_path = data_dict['image_path']
            image_size = data_dict['image_size']
            vps_pred = data_dict['vp']
            
            z_vp = vps_pred[0]
            h_vp = np.array([vps_pred[i] for i in range(1, len(vps_pred))])
            focal_pred, manh_vps, confident = vp2focal(z_vp, h_vp, image_size)  # focal_pred can be -1
            if focal_pred == -1:
                continue
            
            bias = np.array(image_size) / 2
            vps_pred = (np.array(vps_pred) * bias + bias).tolist()

            label = label_dict[image_path]
            focal_gt = label[0]
            vps_gt = label[1]
            
            vps_recall = [nearest(item, vps_pred) for item in vps_gt]
            vps_precision = [nearest(item, vps_gt) for item in vps_pred]
            focal_diff = abs(focal_pred - focal_gt)

            vps_recall_list.append(vps_recall)
            vps_precision_list.append(vps_precision)
            focal_diff_list.append(focal_diff)

        json_out = {'vps_recall': vps_recall_list, 'vps_precision': vps_precision_list, 'focal_diff': focal_diff_list}
        with open(save_name, 'w') as w_op:
            json.dump(json_out, w_op)


if __name__ == '__main__':
    # data_list = ['YUD', 'ScanNet', 'SceneCityUrban3D', 'SUNCG']
    data_list = ['ScanNet_aug', 'SceneCityUrban3D_aug', 'SUNCG_aug']
    dir_path = '/n/fs/vl/xg5/workspace/baseline/Horizon_First_VPdetection/dataset'

    for data_name in data_list:
        print(data_name)
        data_file = '/n/fs/vl/xg5/Datasets/' + data_name + '/label/label.txt'
        with open(data_file, 'r') as op: labels = op.readlines()

        file_name = os.path.join(dir_path, data_name + '/data/data.json')
        save_name = os.path.join(dir_path, data_name + '/data/eval.json')
        run(file_name, save_name, labels)



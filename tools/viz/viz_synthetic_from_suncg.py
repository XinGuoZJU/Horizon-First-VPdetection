import os
import json
import matplotlib.pyplot as plt
import numpy as np


def visualize(line_seg, pred_group, save_name, vp=None):
    # line_seg: n_number x 3, pred_group: list  vp: group_num x dim
    fig = plt.figure()

    if vp is not None:
        axis_list = [1e8, -1e8, 1e8, -1e8]
        for item in vp:
            # x
            if item[0] < axis_list[0]:
                axis_list[0] = item[0]
            if item[0] > axis_list[1]:
                axis_list[1] = item[0]
            # y
            if item[1] < axis_list[2]:
                axis_list[2] = item[1]
            if item[1] > axis_list[3]:
                axis_list[3] = item[1]

        axis_list = [(-1) ** (i + 1) * 1 + int(axis_list[i]) for i in range(4)]
    else:
        axis_list = [-10, 10, -10, 10]
    axis_list = [-5, 5, -5, 5]
    plt.axis(axis_list)

    if vp is not None:
        # draw vp
        for point in vp:
            plt.scatter(point[0], point[1], c='k', s=10, zorder=2)

    color_list = ['y', 'b', 'm', 'r', 'c', 'g', 'w', 'k']
    # draw lines
    #######
    log = np.zeros(8)
    #######

    for i in range(len(line_seg)):
        line = line_seg[i]
        group = int(pred_group[i])
        
        #######
        threshold = 20
        if log[group+1] > threshold:
            continue
        log[group+1] += 1
        #######

        if group == -1:
            color = 'k--'
        else:
            color = color_list[group]
        a, b, c = line
        if b == 0:
            y_r = np.arange(axis_list[0], axis_list[1] + 1, 0.1)
            x_r = -np.ones(len(y_r)) * c / a
        else:
            x_r = np.arange(axis_list[0], axis_list[1] + 1, 0.1)
            y_r = (-c - a * x_r) / b
            idx_low = y_r > axis_list[2]
            idx_up = y_r < axis_list[3]
            idx = idx_low * idx_up
            x_r = x_r[idx]
            y_r = y_r[idx]

        plt.plot(x_r, y_r, color, linewidth=0.5, zorder=1)


    plt.savefig(save_name)
    plt.close()


if __name__ == '__main__':
    org_path = '/n/fs/vl/xg5/workspace/baseline/Horizon-First-VPdetection/dataset/YUD/data/data.json'
    save_path = '/n/fs/vl/xg5/workspace/baseline/Horizon-First-VPdetection/dataset/YUD/viz_line'
    os.makedirs(save_path, exist_ok=True)

    with open(org_path, 'r') as f:
        lines = f.readlines()
        for idx, line in enumerate(lines):
            print(idx)
            item = json.loads(line)
            file_name = item['image_path']
            
            os.makedirs(save_path, exist_ok=True)

            group = np.array(item['group'])
            line_seg = np.array(item['line']).tolist()
            vp = item['vp']

            save_name = os.path.join(save_path, file_name)
            print(save_name)
            visualize(line_seg, group, save_name, vp=vp)



import os


def run(dataset_name, idx, save_op):
    if dataset_name == 'YUD':
        index_file = '/n/fs/vl/xg5/Datasets/YUD/label/index_' + str(idx) + '.txt'
        img_type = 'jpg'
    elif dataset_name == 'ScanNet':
        index_file = '/n/fs/vl/xg5/Datasets/ScanNet/label/index_' + str(idx) + '.txt'
        img_type = 'png'
    elif dataset_name == 'SceneCityUrban3D':
        index_file = '/n/fs/vl/xg5/Datasets/SceneCityUrban3D/label/index_' + str(idx) + '.txt'
        img_type = 'png'
    elif dataset_name == 'SUNCG':
        index_file = '/n/fs/vl/xg5/Datasets/SUNCG/label/index_' + str(idx) + '.txt'
        img_type = 'png'
    elif dataset_name == 'ScanNet_aug':
        index_file = '/n/fs/vl/xg5/Datasets/ScanNet_aug/label/index_' + str(idx) + '.txt'
        img_type = 'jpg'
    elif dataset_name == 'SceneCityUrban3D_aug':
        index_file = '/n/fs/vl/xg5/Datasets/SceneCityUrban3D_aug/label/index_' + str(idx) + '.txt'
        img_type = 'jpg'
    elif dataset_name == 'SUNCG_aug':
        index_file = '/n/fs/vl/xg5/Datasets/SUNCG_aug/label/index_' + str(idx) + '.txt'
        img_type = 'png'
    else:
        raise ValueError('No such dataset!')

    data_path = '/n/fs/vl/xg5/workspace/baseline/Horizon_First_VPdetection/dataset/' \
                                                        + dataset_name + '/output'
    error_file = '/n/fs/vl/xg5/workspace/baseline/Horizon_First_VPdetection/logs/' \
                                    + dataset_name + '_' + str(idx) + '_error.txt'

    file_list = []
    if os.path.isfile(error_file):
        with open(error_file, 'r') as op:
            lines = op.readlines()
            for line in lines:
                line_list = line.strip().split('/')
                image_path = line_list[-2] + '/' + line_list[-1]
                file_list.append(image_path)

    dir_list = os.listdir(data_path)
    for dirs in dir_list:
        dir_path = os.path.join(data_path, dirs)
        sub_dir_list = os.listdir(dir_path)
        for sub_dirs in sub_dir_list:
            sub_dir_path = os.path.join(dir_path, sub_dirs)
            image_path = dirs + '/' + sub_dirs + '.' + img_type
            file_list.append(image_path)

    with open(index_file, 'r') as op:
        lines = op.readlines()
        for line in lines:
            image_name = line.split()[0]
            if image_name not in file_list:
                # print(dataset_name)
                # print(idx)
                # print(image_name)
                # print('\n')
                string = image_name + '\n'
                save_op.write(string)


if __name__ == '__main__':
    # YUD: 1, ScanNet: 265, SceneCityUrban3D: 23, SUNCG: 569
    # data_list = ['YUD', 'ScanNet', 'SceneCityUrban3D', 'SUNCG']
    data_list = ['ScanNet_aug', 'SUNCG_aug', 'SceneCityUrban3D_aug']
    num = 30
    
    for data_name in data_list:
        save_file = 'error_case/' + data_name + '.txt'
        with open(save_file, 'w') as save_op:
            for idx in range(num):
                run(data_name, idx, save_op)

 

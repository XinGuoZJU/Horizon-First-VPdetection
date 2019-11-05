import os


def run(data_name):
    interval_num = 1000
    
    path = 'error_case'
    with open(os.path.join(path, data_name + '.txt'), 'r') as r_op:
        lines = r_op.readlines()
        for i, line in enumerate(lines):
            if i % interval_num == 0:
                ind = int(i / interval_num)
                w_op = open(os.path.join(path, data_name + '_' + str(ind) + '.txt'), 'w')
            w_op.write(line)    
        w_op.close()


if __name__ == '__main__':
    data_list = ['ScanNet', 'SceneCityUrban3D', 'SUNCG']
    for data_name in data_list:
        run(data_name)



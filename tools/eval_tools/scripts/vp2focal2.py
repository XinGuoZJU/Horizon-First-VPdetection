import os
import numpy as np


def vp2focal(z_vp, h_vp, image_size):
    # input should be one zenith vp and multi horizon vps
    # vp: num_vps x 2, [-1,1]
    
    height = image_size[0]
    width = image_size[1]
    # M = -np.matrix([[(height / 2.0) ** 2, 0], [0, (width / 2.0) ** 2]])

    lower_bound = 0.28 * width
    upper_bound = 3.8 * width
    
    infty = 4
    confident = -1
    focal = -1
    manh_vps = []
    accuracy = 2
    nt = 2**(5+min(accuracy, np.log2(width)-5))
    step_h = width / nt
    Dt_h = np.arctan(step_h / width)

    mht_count = 0
    hvp_count =  h_vp.shape[0]
    
    mht = []
    if hvp_count > 1:
        for i in range(hvp_count-1):
            for j in range(i+1, hvp_count):
                vp_i = h_vp[i] * np.array(image_size) / 2
                vp_j = h_vp[j] * np.array(image_size) / 2
                focal_sq = -float((vp_i * vp_j).sum())
                
                if focal_sq > 0:
                    focal = focal_sq ** 0.5
                    
                    if focal > lower_bound and focal < upper_bound:
                        mht_count += 1
                        vv1 = np.append(vp_i, focal)
                        vv2 = np.append(vp_j, focal)
                        vv3 = np.cross(vv1, vv2)
                        vvp_z = np.array([vv3[0] / vv3[2] * focal + height / 2, vv3[1] / vv3[2] * focal + width / 2])

                        at_vz1 = np.arctan(vvp_z[0] / height) / Dt_h
                        at_vz2 = np.arctan((z_vp[0] * height / 2 + height / 2) / height)/ Dt_h
                        score = min(abs(at_vz1 - at_vz2), abs(at_vz1 + at_vz2))

                        mht.append([vvp_z[0], vvp_z[1], focal, i, j, score])

    mht = np.array(mht)
    if mht.shape[0] > 0:
        scores = mht[:, -1]
        index = np.argmin(scores)
        if scores[index] < 2 or abs(z_vp[0] * height / 2) >= infty * height:
            if scores[index] < 2:
                confident = 3
            else:
                index = np.argmax(np.abs(mht[:, 0] - height / 2))
                confident = 1
            
            id1 = int(mht[index, 3])
            id2 = int(mht[index, 4])
            x1 = h_vp[id1][1] * width / 2 + width / 2 
            x2 = h_vp[id2][1] * width / 2 + width / 2
            y = mht[index, 0]
            
            if (y > height / 2 and x1 > x2) or (y < height / 2 and x1 < x2):
                mht[index, 3] = id2
                mht[index, 4] = id1
            focal = mht[index, 2]
            manh_vps = [[h_vp[int(mht[index, 3]), 0], h_vp[int(mht[index, 3]), 1]], 
                        [h_vp[int(mht[index, 4]), 0], h_vp[int(mht[index, 4]), 1]],
                        [z_vp[0], z_vp[1]]]

    if confident == -1 and hvp_count > 0 and abs(z_vp[0] * height / 2) < infty * height:
        id_finite = np.where(np.abs(h_vp[:,1] * width / 2) < infty * width)[0].tolist()
        if len(id_finite) > 0:
            y1 = h_vp[id_finite[0], 0] * height / 2
            y2 = z_vp[0] * height / 2
            focal_sq = float(-y1 * y2)
            if focal_sq > 0:
                focal = focal_sq ** 0.5
                if focal < lower_bound and focal > upper_bound:
                    K = np.matrix([[focal, 0, height / 2], [0, focal, width / 2], [0, 0, 1]])
                    r2 = np.array([h_vp[id_finite[0], 0] * height / 2, h_vp[id_finite[0], 1] * width / 2, focal])
                    r2 = r2 / np.linalg.norm(r2)
                    r3 = np.array([z_vp[0] * height / 2, z_vp[1] * width / 2, focal])
                    r3 = r3 / np.linalg.norm(r3)

                    vpx = K * np.matrix(np.cross(r2, r3)).T
                    vpx = (vpx / vpx[2]).T.tolist()[0]
                    manh_vps = [[vpx[0] / (height / 2), vpx[1] / (width / 2)],
                                [h_vp[id_finite[0], 0], h_vp[id_finite[0], 1]],
                                [z_vp[0], z_vp[1]]]
                    confident = 2

    return focal, manh_vps, confident


if __name__=='__main__':
    image_size = [480, 640]
    h_vp = np.array([[0.0035071074493323805, 2.7071398342638653], [-0.029180096484071664, -1.7742892779112622]])
    z_vp = np.array([208.91234010222863, -0.8571323015478454])
    
    focal, manh_vps, confident = vp2focal(z_vp, h_vp, image_size)
    print(focal)


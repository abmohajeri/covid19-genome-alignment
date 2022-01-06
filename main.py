import numpy as np
import variables
from genome_hash_table import HashTable

from typing import Collection
from numpy.core.fromnumeric import var
import cv2

parts = []
part_num = 0


def nw(part, genome1, genome2, match=variables.Score.match.value, mismatch=variables.Score.mismatch.value,
       gap=variables.Score.gap.value):
       
    global parts
       
    nx = len(genome1)
    ny = len(genome2)
    f = np.zeros((nx + 1, ny + 1))
    f[:, 0] = np.linspace(0, -nx, nx + 1)
    f[0, :] = np.linspace(0, -ny, ny + 1)
    # Pointers to trace through an optimal alignment.
    p = np.zeros((nx + 1, ny + 1))
    p[:, 0] = 3
    p[0, :] = 4
    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if genome1[i] == genome2[j]:
                t[0] = f[i, j] + match
            else:
                t[0] = f[i, j] + mismatch
            t[1] = f[i, j + 1] + gap
            t[2] = f[i + 1, j] + gap
            tmax = np.max(t)
            f[i + 1, j + 1] = tmax
            if t[0] == tmax:
                p[i + 1, j + 1] += 2
            if t[1] == tmax:
                p[i + 1, j + 1] += 3
            if t[2] == tmax:
                p[i + 1, j + 1] += 4
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if p[i, j] in [2, 5, 6, 9]:
            rx.append(genome1[i - 1])
            ry.append(genome2[j - 1])
            i -= 1
            j -= 1
        elif p[i, j] in [3, 5, 7, 9]:
            rx.append(genome1[i - 1])
            ry.append('-')
            i -= 1
        elif p[i, j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(genome2[j - 1])
            j -= 1
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    
    parts.append((part, rx, ry))
       
    score = f[-1, -1]
    return score


def read_genomes(genome_path):
    genome = []
    with open(genome_path) as f:
        while True:
            # Read from file
            c = f.read(1)
            if not c:
                break
            elif c == '\n':
                continue
            else:
                genome.append(c)
    return genome


def create_genome_arrays():
    genome1 = read_genomes(variables.genome1_path)
    genome2 = read_genomes(variables.genome2_path)
    return genome1, genome2


def create_hash_table(genome1):
    ht = HashTable(2 ** variables.num_seq_bits)
    for i in range(len(genome1) - 5):
        substr = genome1[i:i + 6]
        ht.add(substr, i)
    return ht


# extend seed from back and front
def seed_n_extend(genome1, genome2, start1, start2):
    forward_equality = variables.base_seq_len
    backward_equality = 0
    genome1_forward_subseq = genome1[start1:start1 + forward_equality]
    seed_forward = genome2[start2:start2 + forward_equality]
    genome1_backward_subseq = genome1[start1:start1 + forward_equality]
    seed_backward = genome2[start2:start2 + forward_equality]

    while genome1_forward_subseq == seed_forward or genome1_backward_subseq == seed_backward:
        equality = forward_equality + 1
        backward_equality = backward_equality + 1
        genome1_forward_subseq = genome1[start1:start1 + equality]
        seed_forward = genome2[start2:start2 + equality]
        genome1_backward_subseq = genome1[start1 - backward_equality:start1 + equality]
        seed_backward = genome2[start2 - backward_equality:start2 + equality]
        if forward_equality - 1 >= variables.best_cut_threshold:
            return True, start1, start2
        elif backward_equality - 1 >= variables.best_cut_threshold - variables.base_seq_len:
            return True, start1 - backward_equality, start2 - backward_equality
    return False, -1, -1


def reboot_best_cut(genome2, genome2_mid):
    seed = genome2[genome2_mid + 1:genome2_mid + variables.base_seq_len + 1]
    genome2_mid = genome2_mid + 1
    return seed, genome2_mid


def find_best_cut(genome1, genome2):
    genome2_mid = int(len(genome2) / 2)
    initial_genome2_cut = genome2_mid
    initial_genome1_cut = len(genome1) / 2
    seed = genome2[genome2_mid:genome2_mid + variables.base_seq_len]

    while genome2_mid <= len(genome2) - variables.base_seq_len:
        found, genome1_mid = ht.get(seed, len(genome1))
        if found:
            best_cunt_found, genome1_cut, genome2_cut = seed_n_extend(genome1, genome2, genome1_mid, genome2_mid)
            if best_cunt_found:
                return int(genome1_cut), int(genome2_cut)
            else:
                seed, genome2_mid = reboot_best_cut(genome2, genome2_mid)
                continue
        else:
            seed, genome2_mid = reboot_best_cut(genome2, genome2_mid)
            continue
    return int(initial_genome1_cut), int(initial_genome2_cut)


def divide_n_conquer_seq_align(part, genome1, genome2):
    global total_score
    global part_num
       
    p1 = part_num + 1
    part_num += 1
    
    p2 = part_num + 1
    part_num += 1

    if len(genome1) <= variables.minimum_subseq_len and len(genome2) <= variables.minimum_subseq_len:
        return nw(part, genome1, genome2)
    else:
        genome1_cut, genome2_cut = find_best_cut(genome1, genome2)
        return divide_n_conquer_seq_align(p1, genome1[0:genome1_cut], genome2[0:genome2_cut]) + \
                      divide_n_conquer_seq_align(p2, genome1[genome1_cut:-1], genome2[genome2_cut:-1])


genome1, genome2 = create_genome_arrays()
ht = create_hash_table(genome1)


def vconcat_resize(img_list, resize_type='min', interpolation = cv2.INTER_CUBIC):
    w_min = None
    if resize_type == 'min':
        w_min = min(img.shape[1] for img in img_list)
    elif resize_type == 'max':
        w_min = max(img.shape[1] for img in img_list)
    else:
        w_min = max(img.shape[1] for img in img_list)


    im_list_resize = [cv2.resize(img, (w_min, int(img.shape[0] * w_min / img.shape[1])), interpolation = interpolation) for img in img_list]

    return cv2.vconcat(im_list_resize)
  
def putChar(img, char, position, COLORS, font_scale, thickness):
    cv2.rectangle(img, (position[0]-3, position[1]-15), (position[0]+15, position[1]+3), COLORS[char], thickness=-1)
    cv2.rectangle(img, (position[0]-3, position[1]-15), (position[0]+15, position[1]+3),(0,0,0), thickness=1)
    cv2.putText(img, char, position, cv2.FONT_HERSHEY_SIMPLEX, fontScale=font_scale, color=(0, 0, 0), thickness=thickness, lineType =cv2.LINE_AA)
    
def plot_result(parts, gene_parts=variables.gene_parts, MARGIN=variables.MARGIN, 
            CHARACTER_SPACE=variables.CHARACTER_SPACE, font_scale=variables.font_scale, 
            COLORS=variables.COLORS, thickness=variables.thickness, LINE_HEIGHT=variables.LINE_HEIGHT):
    
    sub_plots = []
    parts_sorted = sorted(parts, key = lambda x: x[0])


    for i in range(len(parts_sorted)):

        seqAlign1 = parts_sorted[i][1]
        seqAlign2 = parts_sorted[i][2]

        length = len(seqAlign1)

        count1 = length - seqAlign1.count('-')
        count2 = length - seqAlign2.count('-')


        img = np.zeros([ gene_parts*20 + 2*MARGIN - 15, length*15 + 3*MARGIN,3], dtype=np.uint8)
        img.fill(255)

        for j in range(length):
            
            position1 = (MARGIN + int(j * CHARACTER_SPACE * font_scale), MARGIN)
            position2 = (MARGIN + int(j * CHARACTER_SPACE * font_scale), MARGIN + LINE_HEIGHT)

            putChar(img, seqAlign1[j], position1, COLORS, font_scale, thickness)
            putChar(img, seqAlign2[j], position2, COLORS, font_scale, thickness)

        cv2.putText(img, str(count1), (MARGIN + int(length * CHARACTER_SPACE * font_scale) + 25, MARGIN), cv2.FONT_HERSHEY_SIMPLEX, fontScale=font_scale, color=(0, 0, 0), thickness=thickness, lineType =cv2.LINE_AA)
        cv2.putText(img, str(count2), (MARGIN + int(length * CHARACTER_SPACE * font_scale) + 25, MARGIN + LINE_HEIGHT), cv2.FONT_HERSHEY_SIMPLEX, fontScale=font_scale, color=(0, 0, 0), thickness=thickness, lineType =cv2.LINE_AA)

        sub_plots.append(img)

    result = vconcat_resize(sub_plots, 'max')
    cv2.imwrite('RESULT.jpg', result)




def main():
    print(divide_n_conquer_seq_align(part_num, genome1, genome2))
    plot_result(parts)


if __name__ == '__main__':
    main()

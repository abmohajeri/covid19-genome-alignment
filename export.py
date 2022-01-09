import cv2
import numpy as np
import variables


def vconcat_resize(img_list, resize_type='min', interpolation=cv2.INTER_CUBIC):
    if resize_type == 'min':
        w_min = min(img.shape[1] for img in img_list)
    elif resize_type == 'max':
        w_min = max(img.shape[1] for img in img_list)
    else:
        w_min = max(img.shape[1] for img in img_list)

    im_list_resize = [cv2.resize(img, (w_min, int(img.shape[0] * w_min / img.shape[1])), interpolation=interpolation)
                      for img in img_list]

    return cv2.vconcat(im_list_resize)


def putChar(img, char, position, COLORS, font_scale, thickness):
    cv2.rectangle(img, (position[0] - 3, position[1] - 15), (position[0] + 15, position[1] + 3), COLORS[char],
                  thickness=-1)
    cv2.rectangle(img, (position[0] - 3, position[1] - 15), (position[0] + 15, position[1] + 3), (0, 0, 0), thickness=1)
    cv2.putText(img, char, position, cv2.FONT_HERSHEY_SIMPLEX, fontScale=font_scale, color=(0, 0, 0),
                thickness=thickness, lineType=cv2.LINE_AA)


def plot_result(parts, gene_parts=variables.gene_parts, MARGIN=variables.margin,
                CHARACTER_SPACE=variables.character_space, font_scale=variables.font_scale,
                COLORS=variables.colors, thickness=variables.thickness, LINE_HEIGHT=variables.line_height):
    sub_plots = []
    parts_sorted = sorted(parts, key=lambda x: x[0])

    for i in range(len(parts_sorted)):

        seqAlign1 = parts_sorted[i][1]
        seqAlign2 = parts_sorted[i][2]

        length = len(seqAlign1)

        count1 = length - seqAlign1.count('-')
        count2 = length - seqAlign2.count('-')

        img = np.zeros([gene_parts * 20 + 2 * MARGIN - 15, length * 15 + 3 * MARGIN, 3], dtype=np.uint8)
        img.fill(255)

        for j in range(length):
            position1 = (MARGIN + int(j * CHARACTER_SPACE * font_scale), MARGIN)
            position2 = (MARGIN + int(j * CHARACTER_SPACE * font_scale), MARGIN + LINE_HEIGHT)

            putChar(img, seqAlign1[j], position1, COLORS, font_scale, thickness)
            putChar(img, seqAlign2[j], position2, COLORS, font_scale, thickness)

        cv2.putText(img, str(count1), (MARGIN + int(length * CHARACTER_SPACE * font_scale) + 25, MARGIN),
                    cv2.FONT_HERSHEY_SIMPLEX, fontScale=font_scale, color=(0, 0, 0), thickness=thickness,
                    lineType=cv2.LINE_AA)
        cv2.putText(img, str(count2), (MARGIN + int(length * CHARACTER_SPACE * font_scale) + 25, MARGIN + LINE_HEIGHT),
                    cv2.FONT_HERSHEY_SIMPLEX, fontScale=font_scale, color=(0, 0, 0), thickness=thickness,
                    lineType=cv2.LINE_AA)

        sub_plots.append(img)

    result = vconcat_resize(sub_plots, 'max')
    cv2.imwrite('RESULT.jpg', result)
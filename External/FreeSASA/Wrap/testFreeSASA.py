#  Copyright (c) 2017, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#


import unittest
import os,sys, copy

from rdkit.Chem import rdFreeSASA
from rdkit import Chem

expected = [
    [0, 1, "Polar", 1.64],    [1, 0, "Apolar", 1.88],
    [2, 0, "Apolar", 1.61],   [3, 1, "Polar", 1.42],
    [4, 0, "Apolar", 1.88],   [5, 0, "Apolar", 1.88],
    [6, 1, "Polar", 1.77],    [7, 0, "Apolar", 1.88],
    [8, 1, "Polar", 1.64],    [9, 0, "Apolar", 1.88],
    [10, 0, "Apolar", 1.61],  [11, 1, "Polar", 1.42],
    [12, 0, "Apolar", 1.88],  [13, 0, "Apolar", 1.88],
    [14, 0, "Apolar", 1.61],  [15, 1, "Polar", 1.42],
    [16, 1, "Polar", 1.64],   [17, 1, "Polar", 1.64],
    [18, 0, "Apolar", 1.88],  [19, 0, "Apolar", 1.61],
    [20, 1, "Polar", 1.42],   [21, 0, "Apolar", 1.88],
    [22, 0, "Apolar", 1.88],  [23, 0, "Apolar", 1.88],
    [24, 0, "Apolar", 1.88],  [25, 1, "Polar", 1.64],
    [26, 0, "Apolar", 1.88],  [27, 0, "Apolar", 1.61],
    [28, 1, "Polar", 1.42],   [29, 0, "Apolar", 1.88],
    [30, 0, "Apolar", 1.61],  [31, 0, "Apolar", 1.76],
    [32, 0, "Apolar", 1.76],  [33, 0, "Apolar", 1.76],
    [34, 0, "Apolar", 1.76],  [35, 0, "Apolar", 1.76],
    [36, 1, "Polar", 1.64],   [37, 0, "Apolar", 1.88],
    [38, 0, "Apolar", 1.61],  [39, 1, "Polar", 1.42],
    [40, 0, "Apolar", 1.88],  [41, 0, "Apolar", 1.88],
    [42, 0, "Apolar", 1.88],  [43, 1, "Polar", 1.64],
    [44, 0, "Apolar", 1.88],  [45, 0, "Apolar", 1.61],
    [46, 1, "Polar", 1.42],   [47, 0, "Apolar", 1.88],
    [48, 0, "Apolar", 1.88],  [49, 0, "Apolar", 1.88],
    [50, 0, "Apolar", 1.88],  [51, 1, "Polar", 1.64],
    [52, 1, "Polar", 1.64],   [53, 0, "Apolar", 1.88],
    [54, 0, "Apolar", 1.61],  [55, 1, "Polar", 1.42],
    [56, 0, "Apolar", 1.88],  [57, 1, "Polar", 1.46],
    [58, 0, "Apolar", 1.88],  [59, 1, "Polar", 1.64],
    [60, 0, "Apolar", 1.88],  [61, 0, "Apolar", 1.61],
    [62, 1, "Polar", 1.42],   [63, 0, "Apolar", 1.88],
    [64, 0, "Apolar", 1.88],  [65, 0, "Apolar", 1.88],
    [66, 0, "Apolar", 1.88],  [67, 1, "Polar", 1.64],
    [68, 0, "Apolar", 1.88],  [69, 0, "Apolar", 1.61],
    [70, 1, "Polar", 1.42],   [71, 0, "Apolar", 1.88],
    [72, 1, "Polar", 1.46],   [73, 0, "Apolar", 1.88],
    [74, 1, "Polar", 1.64],   [75, 0, "Apolar", 1.88],
    [76, 0, "Apolar", 1.61],  [77, 1, "Polar", 1.42],
    [78, 1, "Polar", 1.64],   [79, 0, "Apolar", 1.88],
    [80, 0, "Apolar", 1.61],  [81, 1, "Polar", 1.42],
    [82, 0, "Apolar", 1.88],  [83, 0, "Apolar", 1.88],
    [84, 0, "Apolar", 1.88],  [85, 0, "Apolar", 1.88],
    [86, 1, "Polar", 1.64],   [87, 1, "Polar", 1.64],
    [88, 0, "Apolar", 1.88],  [89, 0, "Apolar", 1.61],
    [90, 1, "Polar", 1.42],   [91, 0, "Apolar", 1.88],
    [92, 1, "Polar", 1.46],   [93, 0, "Apolar", 1.88],
    [94, 1, "Polar", 1.64],   [95, 0, "Apolar", 1.88],
    [96, 0, "Apolar", 1.61],  [97, 1, "Polar", 1.42],
    [98, 0, "Apolar", 1.88],  [99, 0, "Apolar", 1.88],
    [100, 0, "Apolar", 1.88], [101, 0, "Apolar", 1.88],
    [102, 1, "Polar", 1.64],  [103, 0, "Apolar", 1.88],
    [104, 0, "Apolar", 1.61], [105, 1, "Polar", 1.42],
    [106, 0, "Apolar", 1.88], [107, 1, "Polar", 1.46],
    [108, 0, "Apolar", 1.88], [109, 1, "Polar", 1.64],
    [110, 0, "Apolar", 1.88], [111, 0, "Apolar", 1.61],
    [112, 1, "Polar", 1.42],  [113, 0, "Apolar", 1.88],
    [114, 0, "Apolar", 1.88], [115, 0, "Apolar", 1.88],
    [116, 0, "Apolar", 1.88], [117, 1, "Polar", 1.64],
    [118, 0, "Apolar", 1.88], [119, 0, "Apolar", 1.61],
    [120, 1, "Polar", 1.42],  [121, 0, "Apolar", 1.88],
    [122, 0, "Apolar", 1.88], [123, 0, "Apolar", 1.61],
    [124, 1, "Polar", 1.42],  [125, 1, "Polar", 1.46],
    [126, 1, "Polar", 1.64],  [127, 0, "Apolar", 1.88],
    [128, 0, "Apolar", 1.61], [129, 1, "Polar", 1.42],
    [130, 0, "Apolar", 1.88], [131, 0, "Apolar", 1.88],
    [132, 0, "Apolar", 1.88], [133, 1, "Polar", 1.64],
    [134, 0, "Apolar", 1.88], [135, 0, "Apolar", 1.61],
    [136, 1, "Polar", 1.42],  [137, 0, "Apolar", 1.88],
    [138, 0, "Apolar", 1.88], [139, 0, "Apolar", 1.61],
    [140, 1, "Polar", 1.42],  [141, 1, "Polar", 1.46],
    [142, 1, "Polar", 1.64],  [143, 0, "Apolar", 1.88],
    [144, 0, "Apolar", 1.61], [145, 1, "Polar", 1.42],
    [146, 0, "Apolar", 1.88], [147, 0, "Apolar", 1.88],
    [148, 0, "Apolar", 1.88], [149, 1, "Polar", 1.64],
    [150, 0, "Apolar", 1.88], [151, 0, "Apolar", 1.61],
    [152, 1, "Polar", 1.42],  [153, 0, "Apolar", 1.88],
    [154, 1, "Polar", 1.46],  [155, 1, "Polar", 1.64],
    [156, 0, "Apolar", 1.88], [157, 0, "Apolar", 1.61],
    [158, 1, "Polar", 1.42],  [159, 0, "Apolar", 1.88],
    [160, 0, "Apolar", 1.61], [161, 1, "Polar", 1.42],
    [162, 1, "Polar", 1.46],  [163, 1, "Polar", 1.64],
    [164, 0, "Apolar", 1.88], [165, 0, "Apolar", 1.61],
    [166, 1, "Polar", 1.42],  [167, 0, "Apolar", 1.88],
    [168, 1, "Polar", 1.46],  [169, 0, "Apolar", 1.88],
    [170, 1, "Polar", 1.64],  [171, 0, "Apolar", 1.88],
    [172, 0, "Apolar", 1.61], [173, 1, "Polar", 1.42],
    [174, 0, "Apolar", 1.88], [175, 0, "Apolar", 1.88],
    [176, 0, "Apolar", 1.88], [177, 0, "Apolar", 1.88],
    [178, 1, "Polar", 1.64],  [179, 0, "Apolar", 1.88],
    [180, 0, "Apolar", 1.61], [181, 1, "Polar", 1.42],
    [182, 0, "Apolar", 1.88], [183, 0, "Apolar", 1.88],
    [184, 0, "Apolar", 1.61], [185, 1, "Polar", 1.42],
    [186, 1, "Polar", 1.46],  [187, 1, "Polar", 1.64],
    [188, 0, "Apolar", 1.88], [189, 0, "Apolar", 1.61],
    [190, 1, "Polar", 1.42],  [191, 0, "Apolar", 1.88],
    [192, 0, "Apolar", 1.61], [193, 1, "Polar", 1.42],
    [194, 1, "Polar", 1.64],  [195, 1, "Polar", 1.64],
    [196, 0, "Apolar", 1.88], [197, 0, "Apolar", 1.61],
    [198, 1, "Polar", 1.42],  [199, 0, "Apolar", 1.88],
    [200, 0, "Apolar", 1.88], [201, 0, "Apolar", 1.88],
    [202, 1, "Polar", 1.64],  [203, 0, "Apolar", 1.88],
    [204, 0, "Apolar", 1.61], [205, 1, "Polar", 1.42],
    [206, 0, "Apolar", 1.88], [207, 0, "Apolar", 1.88],
    [208, 0, "Apolar", 1.88], [209, 0, "Apolar", 1.88],
    [210, 1, "Polar", 1.64],  [211, 1, "Polar", 1.64],
    [212, 0, "Apolar", 1.88], [213, 0, "Apolar", 1.61],
    [214, 1, "Polar", 1.42],  [215, 0, "Apolar", 1.88],
    [216, 1, "Polar", 1.64],  [217, 0, "Apolar", 1.88],
    [218, 0, "Apolar", 1.61], [219, 1, "Polar", 1.42],
    [220, 0, "Apolar", 1.88], [221, 0, "Apolar", 1.88],
    [222, 0, "Apolar", 1.88], [223, 0, "Apolar", 1.88],
    [224, 1, "Polar", 1.64],  [225, 1, "Polar", 1.64],
    [226, 0, "Apolar", 1.88], [227, 0, "Apolar", 1.61],
    [228, 1, "Polar", 1.42],  [229, 0, "Apolar", 1.88],
    [230, 0, "Apolar", 1.88], [231, 0, "Apolar", 1.88],
    [232, 0, "Apolar", 1.88], [233, 1, "Polar", 1.64],
    [234, 0, "Apolar", 1.88], [235, 0, "Apolar", 1.61],
    [236, 1, "Polar", 1.42],  [237, 0, "Apolar", 1.88],
    [238, 0, "Apolar", 1.88], [239, 0, "Apolar", 1.61],
    [240, 1, "Polar", 1.42],  [241, 1, "Polar", 1.64],
    [242, 1, "Polar", 1.64],  [243, 0, "Apolar", 1.88],
    [244, 0, "Apolar", 1.61], [245, 1, "Polar", 1.42],
    [246, 0, "Apolar", 1.88], [247, 0, "Apolar", 1.61],
    [248, 1, "Polar", 1.42],  [249, 1, "Polar", 1.46],
    [250, 1, "Polar", 1.64],  [251, 0, "Apolar", 1.88],
    [252, 0, "Apolar", 1.61], [253, 1, "Polar", 1.42],
    [254, 0, "Apolar", 1.88], [255, 0, "Apolar", 1.88],
    [256, 0, "Apolar", 1.88], [257, 0, "Apolar", 1.88],
    [258, 1, "Polar", 1.64],  [259, 1, "Polar", 1.64],
    [260, 0, "Apolar", 1.88], [261, 0, "Apolar", 1.61],
    [262, 1, "Polar", 1.42],  [263, 0, "Apolar", 1.88],
    [264, 0, "Apolar", 1.88], [265, 0, "Apolar", 1.61],
    [266, 1, "Polar", 1.42],  [267, 1, "Polar", 1.46],
    [268, 1, "Polar", 1.64],  [269, 0, "Apolar", 1.88],
    [270, 0, "Apolar", 1.61], [271, 1, "Polar", 1.42],
    [272, 1, "Polar", 1.64],  [273, 0, "Apolar", 1.88],
    [274, 0, "Apolar", 1.61], [275, 1, "Polar", 1.42],
    [276, 0, "Apolar", 1.88], [277, 0, "Apolar", 1.88],
    [278, 0, "Apolar", 1.88], [279, 0, "Apolar", 1.88],
    [280, 1, "Polar", 1.64],  [281, 0, "Apolar", 1.88],
    [282, 0, "Apolar", 1.61], [283, 1, "Polar", 1.42],
    [284, 0, "Apolar", 1.88], [285, 0, "Apolar", 1.88],
    [286, 0, "Apolar", 1.88], [287, 1, "Polar", 1.64],
    [288, 0, "Apolar", 1.88], [289, 0, "Apolar", 1.61],
    [290, 1, "Polar", 1.42],  [291, 0, "Apolar", 1.88],
    [292, 0, "Apolar", 1.88], [293, 0, "Apolar", 1.88],
    [294, 1, "Polar", 1.64],  [295, 0, "Apolar", 1.88],
    [296, 0, "Apolar", 1.61], [297, 1, "Polar", 1.42],
    [298, 0, "Apolar", 1.88], [299, 0, "Apolar", 1.61],
    [300, 1, "Polar", 1.42],  [301, 1, "Polar", 1.46],
    [302, 1, "Polar", 1.64],  [303, 0, "Apolar", 1.88],
    [304, 0, "Apolar", 1.61], [305, 1, "Polar", 1.42],
    [306, 0, "Apolar", 1.88], [307, 0, "Apolar", 1.88],
    [308, 0, "Apolar", 1.61], [309, 1, "Polar", 1.42],
    [310, 1, "Polar", 1.64],  [311, 1, "Polar", 1.64],
    [312, 0, "Apolar", 1.88], [313, 0, "Apolar", 1.61],
    [314, 1, "Polar", 1.42],  [315, 0, "Apolar", 1.88],
    [316, 0, "Apolar", 1.88], [317, 0, "Apolar", 1.61],
    [318, 1, "Polar", 1.42],  [319, 1, "Polar", 1.64],
    [320, 1, "Polar", 1.64],  [321, 0, "Apolar", 1.88],
    [322, 0, "Apolar", 1.61], [323, 1, "Polar", 1.42],
    [324, 0, "Apolar", 1.88], [325, 0, "Apolar", 1.88],
    [326, 0, "Apolar", 1.88], [327, 1, "Polar", 1.64],
    [328, 0, "Apolar", 1.61], [329, 1, "Polar", 1.64],
    [330, 1, "Polar", 1.64],  [331, 1, "Polar", 1.64],
    [332, 0, "Apolar", 1.88], [333, 0, "Apolar", 1.61],
    [334, 1, "Polar", 1.42],  [335, 0, "Apolar", 1.88],
    [336, 0, "Apolar", 1.88], [337, 0, "Apolar", 1.88],
    [338, 0, "Apolar", 1.88], [339, 1, "Polar", 1.64],
    [340, 0, "Apolar", 1.88], [341, 0, "Apolar", 1.61],
    [342, 1, "Polar", 1.42],  [343, 0, "Apolar", 1.88],
    [344, 0, "Apolar", 1.88], [345, 0, "Apolar", 1.88],
    [346, 0, "Apolar", 1.88], [347, 1, "Polar", 1.64],
    [348, 0, "Apolar", 1.88], [349, 0, "Apolar", 1.61],
    [350, 1, "Polar", 1.42],  [351, 0, "Apolar", 1.88],
    [352, 0, "Apolar", 1.61], [353, 0, "Apolar", 1.76],
    [354, 0, "Apolar", 1.76], [355, 0, "Apolar", 1.76],
    [356, 0, "Apolar", 1.76], [357, 0, "Apolar", 1.76],
    [358, 1, "Polar", 1.64],  [359, 0, "Apolar", 1.88],
    [360, 0, "Apolar", 1.61], [361, 1, "Polar", 1.42],
    [362, 0, "Apolar", 1.88], [363, 1, "Polar", 1.64],
    [364, 0, "Apolar", 1.88], [365, 0, "Apolar", 1.61],
    [366, 1, "Polar", 1.42],  [367, 1, "Polar", 1.64],
    [368, 0, "Apolar", 1.88], [369, 0, "Apolar", 1.61],
    [370, 1, "Polar", 1.42],  [371, 0, "Apolar", 1.88],
    [372, 0, "Apolar", 1.88], [373, 0, "Apolar", 1.88],
    [374, 0, "Apolar", 1.88], [375, 1, "Polar", 1.64],
    [376, 1, "Polar", 1.64],  [377, 0, "Apolar", 1.88],
    [378, 0, "Apolar", 1.61], [379, 1, "Polar", 1.42],
    [380, 0, "Apolar", 1.88], [381, 0, "Apolar", 1.88],
    [382, 0, "Apolar", 1.61], [383, 1, "Polar", 1.42],
    [384, 1, "Polar", 1.64],  [385, 1, "Polar", 1.64],
    [386, 0, "Apolar", 1.88], [387, 0, "Apolar", 1.61],
    [388, 1, "Polar", 1.42],  [389, 0, "Apolar", 1.88],
    [390, 0, "Apolar", 1.88], [391, 0, "Apolar", 1.88],
    [392, 0, "Apolar", 1.88], [393, 1, "Polar", 1.64],
    [394, 0, "Apolar", 1.88], [395, 0, "Apolar", 1.61],
    [396, 1, "Polar", 1.42],  [397, 0, "Apolar", 1.88],
    [398, 0, "Apolar", 1.88], [399, 0, "Apolar", 1.61],
    [400, 1, "Polar", 1.42],  [401, 1, "Polar", 1.46],
    [402, 1, "Polar", 1.64],  [403, 0, "Apolar", 1.88],
    [404, 0, "Apolar", 1.61], [405, 1, "Polar", 1.42],
    [406, 0, "Apolar", 1.88], [407, 0, "Apolar", 1.61],
    [408, 1, "Polar", 1.42],  [409, 1, "Polar", 1.46],
    [410, 1, "Polar", 1.64],  [411, 0, "Apolar", 1.88],
    [412, 0, "Apolar", 1.61], [413, 1, "Polar", 1.42],
    [414, 1, "Polar", 1.64],  [415, 0, "Apolar", 1.88],
    [416, 0, "Apolar", 1.61], [417, 1, "Polar", 1.42],
    [418, 0, "Apolar", 1.88], [419, 0, "Apolar", 1.88],
    [420, 0, "Apolar", 1.88], [421, 1, "Polar", 1.64],
    [422, 0, "Apolar", 1.61], [423, 1, "Polar", 1.64],
    [424, 1, "Polar", 1.64],  [425, 1, "Polar", 1.64],
    [426, 0, "Apolar", 1.88], [427, 0, "Apolar", 1.61],
    [428, 1, "Polar", 1.42],  [429, 0, "Apolar", 1.88],
    [430, 1, "Polar", 1.46],  [431, 0, "Apolar", 1.88],
    [432, 1, "Polar", 1.64],  [433, 0, "Apolar", 1.88],
    [434, 0, "Apolar", 1.61], [435, 1, "Polar", 1.42],
    [436, 0, "Apolar", 1.88], [437, 0, "Apolar", 1.88],
    [438, 0, "Apolar", 1.88], [439, 0, "Apolar", 1.88],
    [440, 1, "Polar", 1.64],  [441, 0, "Apolar", 1.88],
    [442, 0, "Apolar", 1.61], [443, 1, "Polar", 1.42],
    [444, 0, "Apolar", 1.88], [445, 1, "Polar", 1.46],
    [446, 1, "Polar", 1.64],  [447, 0, "Apolar", 1.88],
    [448, 0, "Apolar", 1.61], [449, 1, "Polar", 1.42],
    [450, 0, "Apolar", 1.88], [451, 0, "Apolar", 1.61],
    [452, 1, "Polar", 1.42],  [453, 1, "Polar", 1.46],
    [454, 1, "Polar", 1.64],  [455, 0, "Apolar", 1.88],
    [456, 0, "Apolar", 1.61], [457, 1, "Polar", 1.42],
    [458, 0, "Apolar", 1.88], [459, 0, "Apolar", 1.61],
    [460, 0, "Apolar", 1.76], [461, 0, "Apolar", 1.76],
    [462, 0, "Apolar", 1.76], [463, 0, "Apolar", 1.76],
    [464, 0, "Apolar", 1.61], [465, 1, "Polar", 1.46],
    [466, 1, "Polar", 1.64],  [467, 0, "Apolar", 1.88],
    [468, 0, "Apolar", 1.61], [469, 1, "Polar", 1.42],
    [470, 0, "Apolar", 1.88], [471, 0, "Apolar", 1.61],
    [472, 1, "Polar", 1.42],  [473, 1, "Polar", 1.64],
    [474, 1, "Polar", 1.64],  [475, 0, "Apolar", 1.88],
    [476, 0, "Apolar", 1.61], [477, 1, "Polar", 1.42],
    [478, 0, "Apolar", 1.88], [479, 0, "Apolar", 1.88],
    [480, 0, "Apolar", 1.88], [481, 0, "Apolar", 1.88],
    [482, 1, "Polar", 1.64],  [483, 0, "Apolar", 1.88],
    [484, 0, "Apolar", 1.61], [485, 1, "Polar", 1.42],
    [486, 0, "Apolar", 1.88], [487, 0, "Apolar", 1.88],
    [488, 0, "Apolar", 1.61], [489, 1, "Polar", 1.42],
    [490, 1, "Polar", 1.64],  [491, 1, "Polar", 1.64],
    [492, 0, "Apolar", 1.88], [493, 0, "Apolar", 1.61],
    [494, 1, "Polar", 1.42],  [495, 0, "Apolar", 1.88],
    [496, 0, "Apolar", 1.88], [497, 0, "Apolar", 1.88],
    [498, 0, "Apolar", 1.88], [499, 1, "Polar", 1.64],
    [500, 1, "Polar", 1.64],  [501, 0, "Apolar", 1.88],
    [502, 0, "Apolar", 1.61], [503, 1, "Polar", 1.42],
    [504, 0, "Apolar", 1.88], [505, 0, "Apolar", 1.88],
    [506, 0, "Apolar", 1.61], [507, 1, "Polar", 1.42],
    [508, 1, "Polar", 1.46],  [509, 1, "Polar", 1.64],
    [510, 0, "Apolar", 1.88], [511, 0, "Apolar", 1.61],
    [512, 1, "Polar", 1.42],  [513, 0, "Apolar", 1.88],
    [514, 1, "Polar", 1.46],  [515, 1, "Polar", 1.64],
    [516, 0, "Apolar", 1.88], [517, 0, "Apolar", 1.61],
    [518, 1, "Polar", 1.42],  [519, 0, "Apolar", 1.88],
    [520, 1, "Polar", 1.46],  [521, 0, "Apolar", 1.88],
    [522, 1, "Polar", 1.64],  [523, 0, "Apolar", 1.88],
    [524, 0, "Apolar", 1.61], [525, 1, "Polar", 1.42],
    [526, 0, "Apolar", 1.88], [527, 0, "Apolar", 1.88],
    [528, 0, "Apolar", 1.88], [529, 0, "Apolar", 1.88],
    [530, 1, "Polar", 1.64],  [531, 0, "Apolar", 1.88],
    [532, 0, "Apolar", 1.61], [533, 1, "Polar", 1.42],
    [534, 0, "Apolar", 1.88], [535, 0, "Apolar", 1.61],
    [536, 1, "Polar", 1.64],  [537, 0, "Apolar", 1.76],
    [538, 0, "Apolar", 1.76], [539, 1, "Polar", 1.64],
    [540, 1, "Polar", 1.64],  [541, 0, "Apolar", 1.88],
    [542, 0, "Apolar", 1.61], [543, 1, "Polar", 1.42],
    [544, 0, "Apolar", 1.88], [545, 0, "Apolar", 1.88],
    [546, 0, "Apolar", 1.88], [547, 0, "Apolar", 1.88],
    [548, 1, "Polar", 1.64],  [549, 0, "Apolar", 1.88],
    [550, 0, "Apolar", 1.61], [551, 1, "Polar", 1.42],
    [552, 0, "Apolar", 1.88], [553, 0, "Apolar", 1.88],
    [554, 0, "Apolar", 1.88], [555, 1, "Polar", 1.64],
    [556, 0, "Apolar", 1.88], [557, 0, "Apolar", 1.61],
    [558, 1, "Polar", 1.42],  [559, 0, "Apolar", 1.88],
    [560, 0, "Apolar", 1.88], [561, 0, "Apolar", 1.88],
    [562, 0, "Apolar", 1.88], [563, 1, "Polar", 1.64],
    [564, 0, "Apolar", 1.88], [565, 0, "Apolar", 1.61],
    [566, 1, "Polar", 1.42],  [567, 0, "Apolar", 1.88],
    [568, 0, "Apolar", 1.88], [569, 0, "Apolar", 1.88],
    [570, 1, "Polar", 1.64],  [571, 0, "Apolar", 1.61],
    [572, 1, "Polar", 1.64],  [573, 1, "Polar", 1.64],
    [574, 1, "Polar", 1.64],  [575, 0, "Apolar", 1.88],
    [576, 0, "Apolar", 1.61], [577, 1, "Polar", 1.42],
    [578, 0, "Apolar", 1.88], [579, 0, "Apolar", 1.88],
    [580, 0, "Apolar", 1.88], [581, 0, "Apolar", 1.88],
    [582, 1, "Polar", 1.64],  [583, 0, "Apolar", 1.88],
    [584, 0, "Apolar", 1.61], [585, 1, "Polar", 1.42],
    [586, 0, "Apolar", 1.88], [587, 0, "Apolar", 1.88],
    [588, 0, "Apolar", 1.88], [589, 1, "Polar", 1.64],
    [590, 0, "Apolar", 1.61], [591, 1, "Polar", 1.64],
    [592, 1, "Polar", 1.64],  [593, 1, "Polar", 1.64],
    [594, 0, "Apolar", 1.88], [595, 0, "Apolar", 1.61],
    [596, 1, "Polar", 1.42],  [597, 1, "Polar", 1.64],
    [598, 0, "Apolar", 1.88], [599, 0, "Apolar", 1.61],
    [600, 1, "Polar", 1.42],  [601, 1, "Polar", 1.46]
    ]

class TestCase(unittest.TestCase) :
    def test_basics(self):
        fname = os.path.join(os.environ["RDBASE"],
                             "External", "FreeSASA", "test_data", "1d3z.pdb")
        mol = Chem.MolFromPDBFile(fname)
        radii = rdFreeSASA.classifyAtoms(mol)
        for atom in mol.GetAtoms():
            self.assertEqual( expected[atom.GetIdx()][3], radii[atom.GetIdx()] )
        leeRichards = 5004.79964427
        shrakerupley = 5000.340175

        sasa = rdFreeSASA.CalcSASA(mol, radii=radii)
        self.assertTrue( (sasa-leeRichards) < 1e-5 )

        opts = rdFreeSASA.SASAOpts(rdFreeSASA.ShrakeRupley, rdFreeSASA.Protor)
        sasa = rdFreeSASA.CalcSASA(mol, radii=radii, opts=opts)
        self.assertTrue( (sasa-shrakerupley) < 1e-5 )

        apolar = rdFreeSASA.CalcSASA(mol, radii, query=rdFreeSASA.MakeFreeSasaAPolarAtomQuery(), opts=opts)
        polar = rdFreeSASA.CalcSASA(mol, radii, query=rdFreeSASA.MakeFreeSasaPolarAtomQuery(), opts=opts)

        self.assertTrue( (polar + apolar - 5000.340175) < 1e-5 )

    def test_opts(self):
        fname = os.path.join(os.environ["RDBASE"],
                             "External", "FreeSASA", "test_data", "1d3z.pdb")
        mol = Chem.MolFromPDBFile(fname)
        radii = rdFreeSASA.classifyAtoms(mol)
        for atom in mol.GetAtoms():
            self.assertEqual( expected[atom.GetIdx()][3], radii[atom.GetIdx()] )
        leeRichards = 5004.79964427
        shrakerupley = 5000.340175
        opts = rdFreeSASA.SASAOpts()
        for alg, res in ( (rdFreeSASA.ShrakeRupley, shrakerupley),
                          (rdFreeSASA.LeeRichards, leeRichards)):
            opts.algorithm = alg
            sasa = rdFreeSASA.CalcSASA(mol, radii=radii, opts=opts)
            self.assertTrue( abs(sasa-res) < 1e-5 )
        leeRichards = 5009.93014166
        shrakerupley = 4977.7709106
        opts = rdFreeSASA.SASAOpts()
        opts.probeRadius = 2.0
        for alg, res in ( (rdFreeSASA.ShrakeRupley, shrakerupley),
                          (rdFreeSASA.LeeRichards, leeRichards)):
            opts.algorithm = alg
            sasa = rdFreeSASA.CalcSASA(mol, radii=radii, opts=opts)
            self.assertTrue( abs(sasa-res) < 1e-5 )


if __name__ == '__main__':
  unittest.main()

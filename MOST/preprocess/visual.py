import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import cv2
import os

# 读取数据




def detect_tissue_pixels(
    image_path,
    output_file,
    pixel,
    threshold  # 可调：>0 表示只要有白点就算组织
):
    """
    复刻 MATLAB DBiT 像素识别脚本

    参数：
    - image_path: 二值图路径（黑=背景，白=组织）
    - output_file: 输出坐标文件
    - pixel: 网格大小（50x50）
    - threshold: 判定阈值（0=只要有白点）
    """

    # =========================
    # 读取图像
    # =========================
    
    if not os.path.exists(image_path):
        raise FileNotFoundError(f"Image not found: {image_path}")
    # 如果没指定输出文件 → 默认放到 output_dir
    output_file = Path(output_file)
    if output_file.is_dir() or str(output_file).endswith(("/", "\\")):
        output_file = output_file / "position.txt"
        # 确保输出目录存在
    os.makedirs(output_file.parent, exist_ok=True)
    
    img = cv2.imread(image_path)
    # 转灰度
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # 二值化（和 MATLAB imbinarize 对应）
    _, BW = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    # 转成 0/1
    BW = (BW > 0).astype(np.uint8)

    # =========================
    # 计算网格参数
    # =========================
    pixel_count = 2 * pixel - 1

    numRows, numCols = BW.shape

    pixel_w = numCols / pixel_count
    pixel_h = numRows / pixel_count

    results = []

    # =========================
    #遍历网格
    # =========================
    for i in range(pixel):  # 行

        y = int(round(2 * i * pixel_h))

        for j in range(pixel):  # 列

            x = int(round(2 * j * pixel_w))

            # 防止越界
            y_end = int(round(y + pixel_h))
            x_end = int(round(x + pixel_w))

            region = BW[y:y_end, x:x_end]

            if region.size == 0:
                continue

            # 统计白色像素
            C = np.sum(region)

            # 判定是否有组织
            if C > threshold:
                results.append(f"{j+1}x{i+1}")  # MATLAB 是从1开始

    # =========================
    #输出文件
    # =========================
    with open(output_file, "w") as f:
        f.write(",".join(results))

    print(f"Detected {len(results)} tissue pixels")
    print(f"Saved to {output_file}")
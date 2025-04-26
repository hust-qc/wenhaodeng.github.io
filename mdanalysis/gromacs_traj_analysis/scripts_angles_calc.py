
# 选择三个原子（index 31717, 31718, 32200）
atom1 = u.select_atoms(f"index 31717")  # O5'
atom2 = u.select_atoms(f"index 31718")  # PA
atom3 = u.select_atoms(f"index 32200")  # O1

# 检查是否正确选择了三个原子
if len(atom1) != 1 or len(atom2) != 1 or len(atom3) != 1:
    raise ValueError("One or more atoms not found. Please check the index numbers.")

# 计算角度
angles = []
for ts in u.trajectory:
    # 计算三个原子之间的向量（展平为一维数组）
    vec1 = (atom2.positions - atom1.positions).flatten()
    vec2 = (atom2.positions - atom3.positions).flatten()
    
    # 计算向量模长
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    
    # 检查模长是否为零（避免除以零）
    if norm1 == 0 or norm2 == 0:
        angles.append(np.nan)  # 如果向量无效，记录为NaN
        continue
    
    # 计算点积并标准化
    cos_angle = np.dot(vec1, vec2) / (norm1 * norm2)
    
    # 限制cos_angle在[-1, 1]范围内以避免数值误差
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    
    # 计算角度（弧度）并转换为度
    angle = np.arccos(cos_angle) * (180.0 / np.pi)
    angles.append(angle)

# 转换为numpy数组
angles = np.array(angles)

# 过滤掉NaN值（如果有）并计算平均角度
valid_angles = angles[~np.isnan(angles)]
if len(valid_angles) > 0:
    average_angle = np.mean(valid_angles)
    print(f"Average Angle between atoms: {average_angle:.3f} degrees")
else:
    print("No valid angles calculated.")

# 绘制角度直方图
plt.figure(figsize=(10, 6))
plt.hist(valid_angles, bins=25, color='purple', edgecolor='black', alpha=0.7, label='Angle Distribution')
plt.xlabel('Angle (degrees)', fontsize=16)
plt.ylabel('Frequency', fontsize=16)
# plt.title('Histogram of Angle between Three Atoms', fontsize=16)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)
plt.tick_params(axis='both', labelsize=12)

# 设置X轴和Y轴范围
plt.xlim(120, 180)
plt.ylim(0, 200)

# 保存图像
plt.savefig('angle_histogram_chaina_o5pa_o1.png', dpi=600)
plt.close()

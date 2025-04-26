


# 选择四个原子（index 31714, 31717, 31718, 32200）
atom1 = u.select_atoms(f"index 31714")  # C5'
atom2 = u.select_atoms(f"index 31717")  # O5'
atom3 = u.select_atoms(f"index 31718")  # PA
atom4 = u.select_atoms(f"index 32200")  # O1

# 检查是否正确选择了四个原子
if len(atom1) != 1 or len(atom2) != 1 or len(atom3) != 1 or len(atom4) != 1:
    raise ValueError("One or more atoms not found. Please check the index numbers.")

# 计算二面角
dihedrals = []
for ts in u.trajectory:
    # 获取四个原子的坐标（展平为一维数组）
    p1 = atom1.positions.flatten()
    p2 = atom2.positions.flatten()
    p3 = atom3.positions.flatten()
    p4 = atom4.positions.flatten()
    
    # 计算两个平面的法向量
    # 第一个平面：atom1, atom2, atom3
    vec1 = p2 - p1
    vec2 = p3 - p2
    normal1 = np.cross(vec1, vec2)
    
    # 第二个平面：atom2, atom3, atom4
    vec3 = p3 - p2
    vec4 = p4 - p3
    normal2 = np.cross(vec3, vec4)
    
    # 检查法向量是否有效（模长不为零）
    norm1 = np.linalg.norm(normal1)
    norm2 = np.linalg.norm(normal2)
    if norm1 == 0 or norm2 == 0:
        dihedrals.append(np.nan)  # 如果法向量无效，记录为NaN
        continue
    
    # 计算法向量夹角的余弦
    cos_dihedral = np.dot(normal1, normal2) / (norm1 * norm2)
    
    # 限制cos_dihedral在[-1, 1]范围内以避免数值误差
    cos_dihedral = np.clip(cos_dihedral, -1.0, 1.0)
    
    # 计算二面角（弧度）并转换为度
    dihedral = np.arccos(cos_dihedral) * (180.0 / np.pi)
    
    # 确定二面角的符号（使用三重积）
    triple_product = np.dot(normal1, np.cross(normal2, p3 - p2))
    if triple_product < 0:
        dihedral = -dihedral
    
    dihedrals.append(dihedral)

# 转换为numpy数组
dihedrals = np.array(dihedrals)

# 过滤掉NaN值并计算平均二面角
valid_dihedrals = dihedrals[~np.isnan(dihedrals)]
if len(valid_dihedrals) > 0:
    average_dihedral = np.mean(valid_dihedrals)
    print(f"Average Dihedral Angle between planes: {average_dihedral:.3f} degrees")
else:
    print("No valid dihedral angles calculated.")

# 绘制二面角直方图
plt.figure(figsize=(10, 6))
plt.hist(valid_dihedrals, bins=50, color='blue', edgecolor='black', alpha=0.7, label='Dihedral Angle Distribution')
plt.xlabel('Dihedral Angle (degrees)', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.title('Histogram of Dihedral Angle', fontsize=16)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)
plt.tick_params(axis='both', labelsize=12)

# 设置X轴和Y轴范围
plt.xlim(-180, 180)
plt.ylim(0, 800)

# 保存图像
plt.savefig('dihedral_histogram_c5p_o5p_pa_o1.png', dpi=600)
plt.close()

import numpy as np
import pandas as pd
import plotly.express as px

# --- 1. 定义基础参数 ---
R = 8.314
T = 298.15
RT = R * T

# 溶剂参数: x=水, y=乙醇, z=乙酸乙酯
V_x, d_x, p_x, h_x = 18.1, 15.5, 16.0, 42.3
V_y, d_y, p_y, h_y = 58.5, 15.8, 8.8, 19.4
V_z, d_z, p_z, h_z = 98.5, 15.8,  5.3,  7.2

# --- 2. 新增：根据缩醛度(AD)计算 PVB 溶解度参数的函数 ---
def get_pvb_hsp(AD):
    """根据缩醛度计算 PVB 的 d, p, h 值"""
    poly_d = -0.02112 * AD + 18.82287
    poly_p = -0.09456 * AD + 15.21316
    poly_h = -0.11071 * AD + 20.24478
    return poly_d, poly_p, poly_h

# --- 3. 封装：寻找匹配溶剂的核心函数 ---
def find_solvents_for_pvb(AD, step=0.05, chi_threshold=0.5):
    """
    给定缩醛度 AD，寻找满足 Chi <= threshold 的所有水/醇/酯配比
    """
    # 3.1 获取当前 AD 的聚合物参数
    poly_d, poly_p, poly_h = get_pvb_hsp(AD)
    
    # 3.2 生成网格 (向量化计算)
    x_vals = np.arange(0, 1 + step, step)
    y_vals = np.arange(0, 1 + step, step)
    X_grid, Y_grid = np.meshgrid(x_vals, y_vals)
    
    X_flat = X_grid.flatten()
    Y_flat = Y_grid.flatten()
    
    valid_mask = (X_flat + Y_flat) <= (1.0 + 1e-6)
    x = X_flat[valid_mask]
    y = Y_flat[valid_mask]
    z = 1.0 - x - y
    z[z < 0] = 0  # 修正极小负数
    
    # 3.3 代入公式计算
    V_mix = V_x*x + V_y*y + V_z*z
    
    phi_x = (V_x * x) / V_mix
    phi_y = (V_y * y) / V_mix
    phi_z = (V_z * z) / V_mix
    
    d_mix = d_x*phi_x + d_y*phi_y + d_z*phi_z
    p_mix = p_x*phi_x + p_y*phi_y + p_z*phi_z
    h_mix = h_x*phi_x + h_y*phi_y + h_z*phi_z
    
    # 计算 Hansen 距离平方和 Chi 参数
    distance_sq = 4*(d_mix - poly_d)**2 + (p_mix - poly_p)**2 + (h_mix - poly_h)**2
    chi = 0.6 * (V_mix / RT) * distance_sq
    
    # 3.4 筛选结果
    solution_mask = chi <= chi_threshold
    
    results = pd.DataFrame({
        'x_Water': x[solution_mask].round(3),
        'y_Ethanol': y[solution_mask].round(3),
        'z_EthylAcetate': z[solution_mask].round(3),
        'Chi': chi[solution_mask].round(4)
    })
    
    return results, poly_d, poly_p, poly_h

# ==========================================
# --- 4. 实际使用与测试 (🌟 SCI 级别绘图美化 🌟) ---
# ==========================================

target_AD = 100
 # <-- 你可以随时修改这个值
# 注意：0.0051 会生成大量点，散点重叠能形成类似面域的效果
step_size = 0.0021 # 将步长提取为变量，方便计算面积

# 调用函数
df_solutions, current_d, current_p, current_h = find_solvents_for_pvb(target_AD, step=step_size)

print(f"【当前设定】PVB 缩醛度 AD = {target_AD}")
print(f"【聚合物HSP】d={current_d:.2f}, p={current_p:.2f}, h={current_h:.2f}")
print(f"【计算结果】共找到 {len(df_solutions)} 个满足条件(Chi<=0.5)的配比组合。\n")
# ---------------------------------------------------------
# 📐 新增：溶解窗口面积计算模块
# ---------------------------------------------------------
if len(df_solutions) > 0:
    # 1. 计算溶解点数量
    count_soluble = len(df_solutions)
    
    # 2. 计算单个网格点的代表面积 (dx * dy)
    area_per_point = step_size * step_size
    
    # 3. 计算溶解区在 (x, y) 投影平面上的绝对面积
    area_soluble = count_soluble * area_per_point
    
    # 4. 计算总相图面积 (x + y <= 1 是一个底1高1的直角三角形，面积为 0.5)
    area_total = 0.5
    
    # 5. 计算百分比
    area_ratio = (area_soluble / area_total) * 100
    
    print(f"【面积计算】步长: {step_size}")
    print(f"【面积计算】溶解点数: {count_soluble}")
    print(f"【面积计算】溶解窗口占比: {area_ratio:.2f}% (相对于整个相图)")
    
    # ---------------------------------------------------------
    # 🎨 绘图部分 (将面积添加到标题中)
    # ---------------------------------------------------------
    print("前 5 个可行配比 (x:水, y:乙醇, z:乙酸乙酯):")
    print(df_solutions.head(5))
    
    fig = px.scatter_ternary(
        df_solutions, 
        a="x_Water", b="y_Ethanol", c="z_EthylAcetate",
        color="Chi", 
        color_continuous_scale="Viridis",
        labels={
            "x_Water": "<b>Water</b>",
            "y_Ethanol": "<b>Ethanol</b>",
            "z_EthylAcetate": "<b>Ethyl Acetate</b>",
            "Chi": "<b>Interaction Parameter (<i>&chi;</i>)</b>"
        }
    )
    
    fig.update_traces(marker=dict(size=4, opacity=1.0, line=dict(width=0)))

    
    # 🌟 3. Colorbar (色条) 精细定制 🌟
    fig.update_coloraxes(
        colorbar=dict(
            title=dict(
                text="<b>Chi Parameter (χ)</b>", 
                font=dict(family="Arial", size=18, color="black"),
                side="top"
            ),
            thickness=40,    # 厚度参数
            len=0.75,        # 长度参数
            x=0.80,          # X 位置稍稍靠右，避免挡住右侧的角
            y=0.5,           # 居中对称
            tickfont=dict(family="Arial", size=14, color="black"), 
            outlinewidth=1,  # 给色条加黑色边框
            outlinecolor="black"
        )
    )

    # 🌟 4. 全局字体与三元轴精细化设定 🌟
    # 提取公共坐标轴配置（黑边框、黑刻度、浅色内部辅助网格）
    axis_config = dict(
        title=dict(font=dict(family="Arial", size=20, color="black")),
        tickfont=dict(family="Arial", size=24, color="black"),
        showgrid=True, gridcolor="rgba(211,211,211,0.5)", gridwidth=1, # 浅灰色辅助网格
        linecolor="black", linewidth=2, showline=True,  # 加粗黑色外边框
        ticks="outside", ticklen=6, tickwidth=2, tickcolor="black" # 专业外部刻度
    )

    fig.update_layout(
        font=dict(family="Arial", color="black"), # 全局字体设定
        paper_bgcolor='white', # 外围纯白背景
        plot_bgcolor='white',  # 内部纯白背景
        
        # SCI 图片标题设定 (可选，也可在论文 LaTeX/Word 中直接加 Figure caption)
        title=dict(
             
            font=dict(family="Arial", size=22, color="black"),
            x=0.5, xanchor='center', y=0.95
        ),
        
        # 深度重写三元坐标系样式
        ternary=dict(
            bgcolor="rgba(0,0,0,0)", # 将中心区域强制设为全透明
            aaxis=axis_config,
            baxis=axis_config,
            caxis=axis_config
        ),
        margin=dict(l=40, r=40, b=80, t=80) # 优化留白
    )
    
    fig.show()
    # 推荐保存方式：
    # fig.write_image("SCI_Ternary_Phase_Diagram.png", scale=3)
else:
    print("在此缩醛度下，未找到满足条件的混合溶剂配比。")
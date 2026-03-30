import numpy as np
import plotly.graph_objects as go

# ==========================================
# 1. 基础物理化学参数
# ==========================================
R = 8.314
T = 298.15
RT = R * T

V_x, d_x, p_x, h_x = 18.1, 15.5, 16.0, 42.3
V_y, d_y, p_y, h_y = 58.5, 15.5, 8.8, 19.4
V_z, d_z, p_z, h_z = 98.5, 15.8, 5.3, 7.2

AD_min = 20
AD_max = 100

# ==========================================
# 2. 生成 3D 笛卡尔坐标网格
# ==========================================
grid_res = 120  
X_cart = np.linspace(0, 1, grid_res)
Y_cart = np.linspace(0, np.sqrt(3)/2, int(grid_res * np.sqrt(3)/2))
Z_AD = np.linspace(AD_min, AD_max, grid_res)

X, Y, Z = np.meshgrid(X_cart, Y_cart, Z_AD, indexing='ij')

z_frac = Y / (np.sqrt(3) / 2)
y_frac = X - np.sqrt(3)/3 * Y
x_frac = 1.0 - y_frac - z_frac

# ==========================================
# 3. 核心优化：符号距离(SDF)拼接法
# ==========================================
sd = np.maximum(-x_frac, np.maximum(-y_frac, -z_frac))

x_c = np.clip(x_frac, 0, 1)
y_c = np.clip(y_frac, 0, 1)
z_c = np.clip(z_frac, 0, 1)
sum_c = x_c + y_c + z_c
x_c /= sum_c
y_c /= sum_c
z_c /= sum_c

V_mix = V_x*x_c + V_y*y_c + V_z*z_c
phi_x = (V_x * x_c) / V_mix
phi_y = (V_y * y_c) / V_mix
phi_z = (V_z * z_c) / V_mix

d_mix = d_x*phi_x + d_y*phi_y + d_z*phi_z
p_mix = p_x*phi_x + p_y*phi_y + p_z*phi_z
h_mix = h_x*phi_x + h_y*phi_y + h_z*phi_z

poly_d = -0.02112 * Z + 18.82287
poly_p = -0.09456 * Z + 15.21316
poly_h = -0.11071 * Z + 20.24478

distance_sq = 4*(d_mix - poly_d)**2 + (p_mix - poly_p)**2 + (h_mix - poly_h)**2
chi_base = 0.6 * (V_mix / RT) * distance_sq

# ==========================================
# 4. 图形学魔术：构建平滑墙壁
# ==========================================
barrier = 0.5 + 500.0 * sd
Chi_grid = np.maximum(chi_base, barrier)

# ==========================================
# 5. 可视化渲染 (🌟 恢复物理光照，增强 3D 深度 🌟)
# ==========================================
fig = go.Figure()

fig.add_trace(go.Isosurface(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=Chi_grid.flatten(),
    isomin=0.0,
    isomax=0.5,
    surface_count=2,
    colorscale='Viridis',
    
    # 🌟 核心推荐 1：反转色带 (可选)
    # 将外壳(0.5)变成深紫色。深色背景上的白色高光能将 3D 立体感提升十倍！
    # 如果您必须保留黄色外壳，请将此行设为 False，光照引擎依然能提供阴影。
    reversescale=False, 
    
    caps=dict(x_show=False, y_show=False, z_show=False),
    opacity=0.9, # 稍微调高不透明度，让阴影更扎实
    
    # 🌟 核心推荐 2：高级物理光照引擎 🌟
    lighting=dict(
        ambient=0.3,   # 降低环境光，让没有朝向光源的暗部变黑，产生强烈对比
        diffuse=0.8,   # 极大地增强漫反射，让曲面产生丝滑的明暗渐变
        specular=0.5,  # 开启高光反射，在曲面凸起处形成类似塑料/玻璃的亮斑
        roughness=0.4, # 控制高光范围，0.4 呈现出极具高级感的哑光磨砂质感
        fresnel=0.2    # 边缘光泽，让三维物体边缘更立体
    ),
    name='Solubility Volume',
    
    # Colorbar 保持您要求的完美参数
    colorbar=dict(
        title=dict(
            text="<b>Chi Parameter (χ)</b>", 
            font=dict(family="Arial", size=22, color="black"), 
            side="top"
        ),
        x=0.75,          
        y=0.4,           
        len=0.75,        
        thickness=20,    
        tickfont=dict(family="Arial", size=16, color="black"), 
        outlinewidth=1,  
        outlinecolor="black"
    )
))

# --- 三棱柱边框 ---
pt_water =[0, 0]
pt_ethanol =[1, 0]
pt_ester =[0.5, np.sqrt(3)/2]

lines_x =[pt_water[0], pt_ethanol[0], pt_ester[0], pt_water[0], pt_water[0], pt_ethanol[0], pt_ethanol[0], pt_ester[0], pt_ester[0]]
lines_y =[pt_water[1], pt_ethanol[1], pt_ester[1], pt_water[1], pt_water[1], pt_ethanol[1], pt_ethanol[1], pt_ester[1], pt_ester[1]]
lines_z =[AD_min, AD_min, AD_min, AD_min, AD_max, AD_max, AD_min, AD_min, AD_max]

fig.add_trace(go.Scatter3d(
    x=lines_x, y=lines_y, z=lines_z,
    mode='lines',
    line=dict(color='black', width=4), # 加粗三棱柱边线
    name='Phase Boundaries',
    hoverinfo='skip'
))

fig.add_trace(go.Scatter3d(
    x=[pt_water[0], pt_ethanol[0], pt_ester[0], pt_water[0]],
    y=[pt_water[1], pt_ethanol[1], pt_ester[1], pt_water[1]],
    z=[AD_max, AD_max, AD_max, AD_max],
    mode='lines',
    line=dict(color='black', width=4),
    showlegend=False,
    hoverinfo='skip'
))

# --- 顶点文本标签 ---
for z_pos in[AD_min, AD_max]:
    fig.add_trace(go.Scatter3d(
        x=[pt_water[0]-0.06, pt_ethanol[0]+0.06, pt_ester[0]],
        y=[pt_water[1]-0.06, pt_ethanol[1]-0.06, pt_ester[1]+0.06],
        z=[z_pos, z_pos, z_pos],
        mode='text',
        text=['Water', 'Ethanol', 'Ethyl Acetate'],
        # 🌟 统一为 Arial，采用纯黑色，增大字号 🌟
        textfont=dict(family="Arial", size=18, color='black'),
        showlegend=False,
        hoverinfo='skip'
    ))

# 🌟 全局场景与坐标轴极致净化 🌟
fig.update_layout(
    # 全局字体统一设为 Arial
    font=dict(family="Arial", color="black"),
    
    # 论文图表通常不加标题（或者放在 Caption 里），若需要可取消注释下面这行
    # title=dict(text="3D Solubility Volume for PVB", font=dict(size=20), x=0.5, xanchor='center'),
    
    # 将背景设置为纯白（去除默认的淡蓝色背景）
    paper_bgcolor='white',
    plot_bgcolor='white',
    margin=dict(l=0, r=0, b=0, t=20), # 缩小边缘留白，让图形更大
    
    scene=dict(
        # 彻底隐藏 X 和 Y 轴及网格面
        xaxis=dict(visible=False, showbackground=False),
        yaxis=dict(visible=False, showbackground=False),
        
        # 精细化 Z 轴设定：去除灰面，保留极简黑色线条
        zaxis=dict(
            visible=True,
            showbackground=False,  # 去除灰色背板
            showgrid=False,        # 去除网格线
            zeroline=False,        # 去除零位线
            showline=True,         # 显示 Z 轴主干线
            linecolor='black',     # 坐标轴颜色纯黑
            linewidth=2,           # 轴线粗细
            ticks='outside',       # 刻度朝外
            tickcolor='black',     # 刻度颜色纯黑
            tickwidth=2,
            ticklen=6,
            title=dict(
                text='<b>Acetalization Degree (AD)</b>', 
                font=dict(family="Arial", size=20, color='black') # 极大的 Z 轴标题
            ),
            tickfont=dict(family="Arial", size=16, color='black') # 极大的 Z 轴数字
        ),
        
        # 视角比例锁定
        aspectmode='manual',
        aspectratio=dict(x=1, y=np.sqrt(3)/2, z=1),
        
        # 初始化一个漂亮的相机视角 (可选，确保截图时每次角度一样)
        camera=dict(
            eye=dict(x=1.5, y=1.5, z=0.5)
        )
    )
)

fig.show()
# 如果想要输出高清无损大图用于论文，可以使用下面这行代码 (需要安装 kaleido 库：pip install -U kaleido)
# fig.write_image("SCI_Figure_3D_Phase_Diagram.png", width=1200, height=1000, scale=3)

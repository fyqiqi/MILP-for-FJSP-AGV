import numpy as np
import os
from gurobipy import GRB
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re

from DataRead import getdata
from FJSPAGVMIP import MIPModel

# 设置 Seaborn 风格和字体（设置支持中文字体）
sns.set(style="whitegrid", context="talk", font_scale=1.1)
# 根据系统配置适当设置字体，这里以微软雅黑为例
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Microsoft YaHei']
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# filename = 'FJSSPinstances/1Bilge and Ulusoy/Jobset07.txt'
filename = 'FJSSPinstances/Jobset01.txt'
Data = getdata(filename)
print('data_j', Data['J'], Data['OJ'])
print('DATA_operations_machines', Data['operations_machines'])
print('DATA_operations_machines', Data['operations_times'])

num_operation = []
for i in Data['J']:
    num_operation.append(Data['OJ'][i][-1])
num_operation_max = np.array(num_operation).max()

time_window = np.zeros(shape=(Data['n'], num_operation_max, Data['m']))
for i in range(len(num_operation)):
    for j in range(num_operation[i]):
        mchForJob = Data['operations_machines'][(i + 1, j + 1)]
        for k in mchForJob:
            time_window[i][j][k - 1] = Data['operations_times'][(i + 1, j + 1, k)]
print(time_window)

n = Data['n']
m = Data['m']
J = Data['J']
OJ = Data['OJ']
A = Data['A']
operations_machines = Data['operations_machines']
operations_times = Data['operations_times']
largeM = Data['largeM']

mipmodel = MIPModel(Data)
mipmodel.setParam('OutputFlag', 1)  # 启用输出
mipmodel.optimize()

if mipmodel.status in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.INTERRUPTED]:
    # 采用 pastel 调色板生成淡亮的颜色
    num_jobs = len(J)
    palette = sns.color_palette("pastel", num_jobs)
    # 构造工件到颜色的映射，确保每个工件都有唯一色调（按工件排序）
    job_colors = {job: palette[i] for i, job in enumerate(sorted(J))}

    # 提取机器调度信息
    machine_schedule = {}
    for j in J:
        for i in Data['OJ'][j]:
            for k in Data['operations_machines'][(j, i)]:
                if mipmodel.getVarByName(f"x_{j}_{i}_{k}").X > 0.5:
                    start = mipmodel.getVarByName(f"s_{j}_{i}").X
                    duration = {key: value * 2 for key, value in Data['operations_times'].items()}
                    machine_schedule.setdefault(k, []).append({
                        'job': j,
                        'operation': i,
                        'start': start,
                        'end': start + duration,
                        'duration': duration
                    })

    # 提取 AGV 运输调度信息
    agv_schedule = {}
    for j in J:
        for i in Data['OJ'][j]:
            if mipmodel.getVarByName(f"eta_{j}_{i}").X > 0.5:
                for r in Data['A']:
                    if mipmodel.getVarByName(f"mu_{j}_{i}_{r}").X > 0.5:
                        start = mipmodel.getVarByName(f"ST_{j}_{i}").X
                        duration = mipmodel.getVarByName(f"Tt_{j}_{i}").X
                        agv_schedule.setdefault(r, []).append({
                            'job': j,
                            'operation': i,
                            'start': start,
                            'end': start + duration,
                            'duration': duration
                        })

    # 绘制甘特图
    plt.figure(figsize=(14, 9))
    y_ticks = []
    y_labels = []

    # 绘制机器调度部分
    for idx, machine in enumerate(sorted(machine_schedule.keys())):
        y_pos = idx * 2
        y_ticks.append(y_pos)
        y_labels.append(f"机器 {machine}")
        for op in machine_schedule[machine]:
            plt.barh(
                y=y_pos,
                width=op['duration'],
                left=op['start'],
                color=job_colors[op['job']],
                edgecolor='gray',
                linewidth=1,
                label=f"Job {op['job']}"
            )

    # 绘制 AGV 调度部分，并稍加透明调整
    agv_base = len(machine_schedule) * 2
    for idx, (agv, tasks) in enumerate(sorted(agv_schedule.items())):
        y_pos = agv_base + idx * 2
        y_ticks.append(y_pos)
        y_labels.append(f"AGV {agv}")
        for task in tasks:
            plt.barh(
                y=y_pos,
                width=task['duration'],
                left=task['start'],
                color=job_colors[task['job']],
                edgecolor='gray',
                linewidth=1,
                alpha=0.8,
                label=f"Job {task['job']}"
            )

    plt.yticks(y_ticks, y_labels,size = 16)

    # 提取所有图例项并按 Job 编号排序（只处理以 "Job" 开头的图例）
    handles, labels = plt.gca().get_legend_handles_labels()
    job_items = []
    pattern = re.compile(r"Job\s*(\d+)")
    for h, l in zip(handles, labels):
        m = pattern.search(l)
        if m:
            job_items.append((int(m.group(1)), l, h))

    # 去重并按编号排序，确保顺序为 Job 1, Job 2, …
    seen = {}
    for num, l, h in job_items:
        if l not in seen:
            seen[l] = (num, h)
    sorted_items = sorted(seen.items(), key=lambda x: x[1][0])

    # 生成排序后的句柄和标签列表
    sorted_labels = [item[0] for item in sorted_items]
    sorted_handles = [item[1][1] for item in sorted_items]

    # 将图例放在图表上方，横向排列，每行显示 4 个图例项
    ncol = 4
    plt.legend(sorted_handles, sorted_labels,
               loc='upper center',
               bbox_to_anchor=(0.5, 1.18),
               ncol=ncol,
               fancybox=True,
               shadow=False)

    plt.xlabel("时间",size=16)
    # plt.title("Gantt Chart of Machine & AGV Scheduling")
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.show()
else:
    print("No solution found")
# # 初始化工件的就绪时间
# job_ready_times = {job: 0 for job in J}
# # 按动态就绪时间生成加工顺序链
# all_operations = []
# for machine, tasks in machine_schedule.items():
#     for task in tasks:
#         all_operations.append({
#             'job': task['job'],
#             'operation': task['operation'],
#             'machine': machine,
#             'start': task['start'],
#             'end': task['end'],
#             'duration': task['duration']
#         })
#
# for agv, tasks in agv_schedule.items():
#     for task in tasks:
#         for op in all_operations:
#             if op['job'] == task['job'] and op['operation'] == task['operation']:
#                 op['agv'] = agv
#                 break

# # 动态生成加工顺序链
# operation_chain = []
# encoded_values = []
# action_set = []
#
# # 初始化上一动作的结束时间
# last_end_time = 0
# # 初始化机器的就绪时间
# machine_ready_times = {machine: 0 for machine in machine_schedule.keys()}
#
# # 重置就绪时间
# job_ready_times = {job: 0 for job in J}
# # 拷贝操作集合用于动态调度
# remaining_ops = all_operations.copy()
#
# while remaining_ops:
#     # 找到具有最小 max(工件就绪时间, 机器就绪时间) 的工序
#     remaining_ops.sort(
#         key=lambda op: ( machine_ready_times[op['machine']], op['start']))
#     # for op in remaining_ops:
#     #     print(op)
#     current_op = remaining_ops.pop(0)
#
#     job_id = current_op['job']
#     op_id = current_op['operation']
#     machine_id = current_op['machine']
#     agv_id = int(current_op.get('agv', 1))  # 默认 AGV ID 为 1
#
#     # 当前工序的实际开始时间为 max(工件就绪时间，机器就绪时间)
#     real_start_time = max(job_ready_times[job_id], machine_ready_times[machine_id])
#     real_end_time = real_start_time + current_op['end']
#
#     # 更新就绪时间
#     machine_ready_times[machine_id] = real_end_time
#
#     # 计算动作编码
#     actions = ((agv_id - 1) * Data['m'] * Data['n']) + ((job_id - 1) * Data['m']) + machine_id - 1
#     chosen_agv = actions // (Data['m'] * Data['n'])
#     chosen_job = (actions - Data['n'] * Data['m'] * chosen_agv) // Data['m']
#     chosen_mch = (actions - Data['n'] * Data['m'] * chosen_agv) % Data['m']
#
#     # 添加到结果链中
#     operation_chain.append(f"O{job_id}{op_id}")
#     encoded_values.append((chosen_agv, chosen_job, chosen_mch))
#     action_set.append(actions)
#
# print("\n加工顺序链：")
# print(" → ".join(operation_chain))
#
# print("\n编码值序列 (AGV, Job, Machine)：")
# print(" → ".join(str(v) for v in encoded_values))
# print(" → ".join(str(v) for v in action_set))

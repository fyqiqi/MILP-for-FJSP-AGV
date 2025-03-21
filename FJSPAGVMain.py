import numpy as np
import time
import os
from map import get_travel_time
from gurobipy import GRB
from matplotlib import pyplot as plt

from DataRead import getdata
from FJSPAGVMIP import MIPModel

filename = './FJSSPinstances/MK/MK05.fjs'
Data = getdata(filename)
print('data_j', Data['J'], Data['OJ'])
print('DATA_operations_machines', Data['operations_machines'])
print('DATA_operations_machines', Data['operations_times'])

num_operation = []
for i in Data['J']:
    num_operation.append(Data['OJ'][i][-1])
print(num_operation)
num_operation_max = np.array(num_operation).max()

time_window = np.zeros(shape=(Data['n'], num_operation_max, Data['m']))

for i in range(len(num_operation)):
    for j in range(num_operation[i]):
        mchForJob = Data['operations_machines'][(i + 1, j + 1)]
        for k in mchForJob:
            time_window[i][j][k - 1] = Data['operations_times'][(i + 1, j + 1, k)]
print(time_window)
# map_window = get_travel_time()
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

# 确保即使没有最优解也能绘制甘特图
if mipmodel.status in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.INTERRUPTED]:
    # 提取机器调度信息
    machine_schedule = {}
    for j in Data['J']:
        for i in Data['OJ'][j]:
            for k in Data['operations_machines'][(j, i)]:
                if mipmodel.getVarByName(f"x_{j}_{i}_{k}").X > 0.5:
                    start = mipmodel.getVarByName(f"s_{j}_{i}").X
                    duration = Data['operations_times'][(j, i, k)]
                    machine_schedule.setdefault(k, []).append({
                        'job': j,
                        'operation': i,
                        'start': start,
                        'end': start + duration,
                        'duration': duration
                    })

    # 提取AGV运输信息
    agv_schedule = {}
    for j in Data['J']:
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
    plt.figure(figsize=(12, 8))

    # 绘制机器调度
    colors = plt.cm.tab20.colors
    y_ticks = []
    y_labels = []

    # 绘制机器部分
    for idx, machine in enumerate(sorted(machine_schedule.keys())):
        y_pos = idx * 2
        y_ticks.append(y_pos)
        y_labels.append(f"Machine {machine}")
        for op in machine_schedule[machine]:
            plt.barh(y=y_pos,
                     width=op['duration'],
                     left=op['start'],
                     color=colors[op['job'] % 20],
                     edgecolor='black',
                     label=f'Job {op["job"]}')

    # 绘制AGV部分
    agv_base = len(machine_schedule) * 2
    for idx, (agv, tasks) in enumerate(agv_schedule.items()):
        y_pos = agv_base + idx * 2
        y_ticks.append(y_pos)
        y_labels.append(f"AGV {agv}")
        for task in tasks:
            plt.barh(y=y_pos,
                     width=task['duration'],
                     left=task['start'],
                     color=colors[task['job'] % 20],
                     edgecolor='black',
                     alpha=0.7)

    plt.xlabel('Time')
    plt.ylabel('Resources')
    plt.yticks(y_ticks, y_labels)
    plt.title('Schedule Gantt Chart')

    # 处理图例去重
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(),
               bbox_to_anchor=(1.05, 1),
               loc='upper left')

    # plt.tight_layout()
    plt.grid(axis='x')
    plt.show()
else:
    print("No solution found")

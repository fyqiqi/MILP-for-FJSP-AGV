from gurobipy import Model, GRB, quicksum
import sys
from map import get_travel_time, read_and_convert_txt


def MIPModel(Data):
    n = Data['n']
    m = Data['m']
    J = Data['J']
    A = Data['A']
    OJ = Data['OJ']
    operations_machines = Data['operations_machines']
    # operations_times = Data['operations_times']*2
    operations_times = {key: value * 3 for key, value in Data['operations_times'].items()}
    largeM = Data['largeM']
    # TT = read_and_convert_txt('./FJSSPinstances/Layout1.txt')
    TT = get_travel_time()/2.0

    model = Model("FJSP_PPF")
    model.Params.OutputFlag = 0
    # Decision Variables
    x = {}      # φ[j,i,k]
    eta = {}    # η[j,i]
    mu = {}     # μ[j,i,r]
    z = {}      # z[j1,i1,j2,i2,r]
    s = {}      # SO[j,i]
    ST = {}     # ST[j,i]
    Tt = {}     # Tt[j,i]
    seq = {}
    msq = {}
    cmax = model.addVar(vtype=GRB.INTEGER, name="cmax")

    for j in J:
        for i in OJ[j]:
            s[j, i] = model.addVar(vtype=GRB.INTEGER, name=f"s_{j}_{i}")
            ST[j, i] = model.addVar(vtype=GRB.INTEGER, name=f"ST_{j}_{i}")
            Tt[j, i] = model.addVar(vtype=GRB.INTEGER, name=f"Tt_{j}_{i}")
            eta[j, i] = model.addVar(vtype=GRB.BINARY, name=f"eta_{j}_{i}")
            for k in operations_machines[j, i]:
                x[j, i, k] = model.addVar(vtype=GRB.BINARY, name=f"x_{j}_{i}_{k}")
            for r in A:
                mu[j, i, r] = model.addVar(vtype=GRB.BINARY, name=f"mu_{j}_{i}_{r}")
    for r in A:
        for j1 in J:
            for j2 in J:
                for i1 in OJ[j1]:
                    for i2 in OJ[j2]:
                        z[j1, i1, j2, i2, r] = model.addVar(vtype=GRB.BINARY, name=f"z_{j1}_{i1}_{j2}_{i2}_{r}")
    for k in range(1, m + 1):
        for j1 in J:
            for j2 in J:
                for i1 in OJ[j1]:
                    for i2 in OJ[j2]:
                        msq[j1, i1, j2, i2, k] = model.addVar(vtype=GRB.BINARY, name=f"z_{j1}_{i1}_{j2}_{i2}_{k}")
    # Objective
    model.setObjective(cmax, GRB.MINIMIZE)

    # Constraints
    # 1. 机器分配 (Sub-problem 1)
    for j in J:
        for i in OJ[j]:
            model.addConstr(
                quicksum(x[j, i, k] for k in operations_machines[j, i]) == 1,
                f"machine_assignment_{j}_{i}"
            )
        # 7. 同一机器上的工序顺序约束（新增）
    for k in range(1, m + 1):  # 遍历所有机器
        # 收集所有可能分配到机器k的工序
        ops_on_k = [(j, i) for (j, i) in operations_machines if k in operations_machines[(j, i)]]
        # 生成所有可能的工序对
        for idx1 in range(len(ops_on_k)):
            for idx2 in range(idx1 + 1, len(ops_on_k)):
                j1, i1 = ops_on_k[idx1]
                j2, i2 = ops_on_k[idx2]
                # 创建顺序变量 sq
                sq = model.addVar(vtype=GRB.BINARY, name=f"sq_{j1}_{i1}_{j2}_{i2}_{k}")

                # 约束1: j1在j2之前
                model.addConstr(
                    s[j2, i2] >= s[j1, i1] + x[j1, i1, k] * operations_times[j1, i1, k]
                    - largeM * (2 - x[j1, i1, k] - x[j2, i2, k]) - largeM * (1 - sq),
                    name=f"machine_order_forward_{j1}_{i1}_{j2}_{i2}_{k}"
                )
                # 约束2: j2在j1之前
                model.addConstr(
                    s[j1, i1] >= s[j2, i2] + x[j2, i2, k] * operations_times[j2, i2, k]
                    - largeM * (2 - x[j1, i1, k] - x[j2, i2, k]) - largeM * sq,
                    name=f"machine_order_backward_{j1}_{i1}_{j2}_{i2}_{k}"
                )
    #todo eta[j,i]的值定义有误
    for j in J:
        # 首工序必须存在运输任务
        first_op = OJ[j][0]
        model.addConstr(eta[j, first_op] == 1, f"eta_first_{j}")
        # 非首工序
        for idx in range(1, len(OJ[j])):
            i = OJ[j][idx]
            prev_i = OJ[j][idx-1]
            model.addConstr(
                eta[j, i] == 1 - quicksum(
                    x[j,i,k]*x[j,prev_i,k]
                    for k in operations_machines[j, prev_i]
                    if k in operations_machines[j, i]
                ),
                f"eta_{j}_{i}"
            )
    #不运输时保证ST与Tt为0
    #不运输不用考虑ST的值，随便给一个
    for j in J:
        for i in OJ[j]:
            model.addConstr(
                (eta[j, i] == 0) >> (ST[j, i] == 0),
                f"ST_zero_when_no_transport_{j}_{i}"
            )
            model.addConstr(
                (eta[j, i] == 0) >> (Tt[j, i] == 0),
                f"Tt_zero_when_no_transport_{j}_{i}"
            )
    # 3. 运输时间计算 (Tt) y为i 是否在k_pre和k_curr上面加工
    for j in J:
        for idx, i in enumerate(OJ[j]):
            if idx == 0:
                model.addConstr(Tt[j, i] == quicksum(x[j, i, k] * TT[0][k] for k in operations_machines[j, i]))
            else:
                prev_i = OJ[j][idx - 1]
                y = {}
                for k_prev in operations_machines[j, prev_i]:
                    for k_curr in operations_machines[j, i]:
                        y[j, i, k_prev, k_curr] = model.addVar(vtype=GRB.BINARY)
                        model.addConstr(y[j, i, k_prev, k_curr] <= x[j, prev_i, k_prev])
                        model.addConstr(y[j, i, k_prev, k_curr] <= x[j, i, k_curr])
                        model.addConstr(y[j, i, k_prev, k_curr] >= x[j, prev_i, k_prev] + x[j, i, k_curr] - 1)
                model.addConstr(
                    Tt[j, i] == quicksum(
                        y[j, i, k_prev, k_curr] * TT[k_prev][k_curr]
                        # (x[j,i,k_curr]*x[j,prev_i,k_prev]) * TT[k_prev][k_curr]
                        for k_prev in operations_machines[j, prev_i]
                        for k_curr in operations_machines[j, i]
                    ),
                    name=f"Tt_{j}_{i}"
                )

    # 4. AGV分配 (Sub-problem 3)
    for j in J:
        for i in OJ[j]:
            model.addConstr(
                quicksum(mu[j, i, r] for r in A) == eta[j, i],
                f"AGV_assignment_{j}_{i}"
            )
    # 定义辅助变量 当且仅当 mu[j1,i1,r] == 1 且 mu[j2,i2,r] == 1 时，w = 1：
    w = {}
    for r in A:
        for j1 in J:
            for j2 in J:
                for i1 in OJ[j1]:
                    for i2 in OJ[j2]:
                        if (j1, i1) != (j2, i2):
                            w[j1, i1, j2, i2, r] = model.addVar(
                                vtype=GRB.BINARY,
                                name=f"w_{j1}_{i1}_{j2}_{i2}_{r}"
                            )
    for r in A:
        for j1 in J:
            for j2 in J:
                for i1 in OJ[j1]:
                    for i2 in OJ[j2]:
                        if  (j1, i1) != (j2, i2):
                            # w <= mu[j1,i1,r]
                            model.addConstr(
                                w[j1, i1, j2, i2, r] <= mu[j1, i1, r],
                                name=f"w_upper1_{j1}_{i1}_{j2}_{i2}_{r}"
                            )
                            # w <= mu[j2,i2,r]
                            model.addConstr(
                                w[j1, i1, j2, i2, r] <= mu[j2, i2, r],
                                name=f"w_upper2_{j1}_{i1}_{j2}_{i2}_{r}"
                            )
                            # w >= mu[j1,i1,r] + mu[j2,i2,r] - 1
                            model.addConstr(
                                w[j1, i1, j2, i2, r] >= mu[j1, i1, r] + mu[j2, i2, r] - 1,
                                name=f"w_lower_{j1}_{i1}_{j2}_{i2}_{r}"
                            )
    # 5. AGV任务排序 (Sub-problem 4)
    #todo  z决策变量的值定义有误
    # 5. AGV任务排序约束（修正后的运输时间计算）
    for r in A:
        for j1 in J:
            for j2 in J:
                for i1 in OJ[j1]:
                    for i2 in OJ[j2]:
                        if (j1, i1) == (j2, i2):
                            continue  # 跳过相同任务
                        # 获取j2的前序工序索引
                        if i2 == OJ[j2][0]:  # 首工序
                            # 首工序的运输起点是depot（0）
                            for k1 in operations_machines[j1, i1]:
                                # AGV从j1,i1的机器k1移动到depot（0），再运输到当前机器k2
                                for k2_curr in operations_machines[j2, i2]:
                                    # 运输时间TT[k1][0] + 当前工序运输时间Tt[j2,i2]（已包含在ST）
                                    model.addConstr(
                                        (z[j1, i1, j2, i2, r] == 1) >> (
                                                ST[j2, i2] >= ST[j1, i1] + Tt[j1, i1] + TT[k1][0]
                                                - largeM * (2 - x[j1, i1, k1] - x[j2, i2, k2_curr])
                                        ),
                                        f"AGV_order_{j1}_{i1}_{j2}_{i2}_first_{r}"
                                    )
                        else:
                            # 非首工序，前序工序的机器k_prev
                            prev_i2 = OJ[j2][OJ[j2].index(i2) - 1]
                            # 遍历可能的前序机器k_prev和当前机器k_curr
                            for k1 in operations_machines[j1, i1]:
                                for k_prev in operations_machines[j2, prev_i2]:
                                    for k_curr in operations_machines[j2, i2]:
                                        # AGV从j1,i1的机器k1移动到前序机器k_prev，再运输到k_curr
                                        model.addConstr(
                                            (z[j1, i1, j2, i2, r] == 1) >> (
                                                    ST[j2, i2] >= ST[j1, i1] + Tt[j1, i1] + TT[k1][k_prev]
                                                    - largeM * (3 - x[j1, i1, k1] - x[j2, prev_i2, k_prev] - x[
                                                j2, i2, k_curr])
                                            ),
                                            f"AGV_order_{j1}_{i1}_{j2}_{i2}_{r}"
                                        )
                        # 确保顺序互斥
                        if (j1, i1) < (j2, i2):  # 避免重复处理
                            model.addConstr(
                                z[j1, i1, j2, i2, r] + z[j2, i2, j1, i1, r] <= 1,
                                name=f"AGV_order_mutex_{j1}_{i1}_{j2}_{i2}_{r}"
                            )
                        if (j1, i1) != (j2, i2):
                            model.addConstr(
                                (w[j1,i1,j2,i2,r] == 1 ) >> (
                                        z[j1, i1, j2, i2, r] + z[j2, i2, j1, i1, r] == 1
                                ),
                                name=f"AGV_order_mutex_{j1}_{i1}_{j2}_{i2}_{r}"
                            )
                            # 条件：至少一个任务未被分配到AGV r
                            model.addConstr(
                                z[j1, i1, j2, i2, r] <= mu[j1, i1, r],
                                name=f"z_upper1_{j1}_{i1}_{j2}_{i2}_{r}"
                            )
                            model.addConstr(
                                z[j1, i1, j2, i2, r] <= mu[j2, i2, r],
                                name=f"z_upper2_{j1}_{i1}_{j2}_{i2}_{r}"
                            )
                            model.addConstr(
                                z[j2, i2, j1, i1, r] <= mu[j1, i1, r],
                                name=f"z_upper1_{j2}_{i2}_{j1}_{i1}_{r}"
                            )
                            model.addConstr(
                                z[j2, i2, j1, i1, r] <= mu[j2, i2, r],
                                name=f"z_upper2_{j2}_{i2}_{j1}_{i1}_{r}"
                            )
    #保证加工顺序
    for j in J:
        for idx in range(1, len(OJ[j])):  # 从第二道工序开始
            i = OJ[j][idx]
            prev_i = OJ[j][idx - 1]
            model.addConstr(
                s[j, i] >= s[j, prev_i] + quicksum(
                    x[j, prev_i, k] * operations_times[j, prev_i, k]
                    for k in operations_machines[j, prev_i]
                ),
                name=f"op_sequence_{j}_{i}"
            )
    #运输时间在加工完成之后
    for j in J:
        for idx in range(1, len(OJ[j])):
            i = OJ[j][idx]
            prev_i = OJ[j][idx - 1]
            model.addConstr(
                ST[j, i] >= s[j, prev_i] + quicksum(
                    x[j, prev_i, k] * operations_times[j, prev_i, k]
                    for k in operations_machines[j, prev_i]
                ) - largeM * (1 - eta[j, i]),
                name=f"ST_after_prev_op_{j}_{i}"
            )

    # 6. 处理与运输时序耦合
    for j in J:
        for idx, i in enumerate(OJ[j]):
            if idx == 0:  # 首工序
                model.addConstr(
                    s[j, i] >= ST[j, i] + Tt[j, i],
                    f"coupling_first_{j}_{i}"
                )
            else:         # 非首工序
                prev_i = OJ[j][idx-1]
                # 当η=1时，处理任务在运输后开始
                model.addConstr(
                    s[j, i] >= ST[j, i] + Tt[j, i] - largeM * (1 - eta[j, i]),
                    f"coupling_transport_{j}_{i}"
                )
                # η=0时：紧接前序工序结束
                model.addConstr(
                    (eta[j, i] == 0) >> (
                            s[j, i] >= s[j, prev_i] + quicksum(
                        x[j, prev_i, k] * operations_times[j, prev_i, k]
                        for k in operations_machines[j, prev_i]
                    )
                    ),
                    f"coupling_no_transport_{j}_{i}"
                )

    # 7. Makespan定义
    for j in J:
        for i in OJ[j]:
            model.addConstr(
                cmax >= s[j, i] + quicksum(
                    x[j, i, k] * operations_times[j, i, k]
                    for k in operations_machines[j, i]
                ),
                f"cmax_{j}"
            )
    model.Params.MIPFocus = 2  # 让求解器专注于尽可能多地改进当前的最优解

    model.Params.TimeLimit = 200000
    model.update()
    return model
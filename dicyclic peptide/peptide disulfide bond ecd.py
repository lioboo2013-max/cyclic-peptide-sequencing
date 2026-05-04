# #这段 Python 代码构建了一个**生物质谱 ECD 碎裂模拟器**。它的核心目的是解决一个复杂的排列组合问题：当多肽含有多个二硫键时，ECD（电子捕获解离）可能打断骨架，也可能打断二硫键，或者两者都打断。

# ### 1\. 基础配置与物理常数模块
#
# 代码开头定义了地基，即“原子有多重”。
#
#   * **`AA_MASS` 字典**: 存储 20 种氨基酸残基的单一同位素质量（Monoisotopic Mass）。例如 `C` (Cysteine) 是 103.00919 Da，这是还原状态（-SH）下的残基质量。
#   * **`CONSTANTS` 字典**:
#       * `PROTON` (1.00728): 质谱观测的是带电离子，所以计算时总要加一个质子质量 ($[M+H]^+$)。
#       * `NH3`/`OH`/`NH`: 用于计算 c 离子和 z 离子的末端基团。
#       * `SS_LOSS` (2 \* 1.00783): **这是关键。** 当两个 Cys 形成一个二硫键时，会脱去 2 个氢原子。代码会在计算质量时减去这个值。
#
# -----
#
# ### 2\. 拓扑结构生成模块 (排列组合)
#
# 这一部分负责回答：“这 6 个半胱氨酸（Cys）到底是怎么连的？”
#
#   * **`get_cys_indices(seq)`**:
#       * **功能**: 扫描序列，找出所有 'C' 的位置（例如：第4、10、14位...）。
#   * **`generate_disulfide_patterns(cys_list)`**:
#       * **逻辑**: 这是一个**递归函数**。
#       * 它采用“固定第一个，遍历剩余伙伴”的策略。比如有 [C1, C2, C3, C4]。
#         1.  固定 C1，它可以连 C2。剩余 [C3, C4] 递归去配对。
#         2.  固定 C1，它可以连 C3。剩余 [C2, C4] 递归去配对。
#         3.  ...以此类推。
#       * **结果**: 对于 6 个 Cys，它会精准生成 15 种连接模式。
#
# -----
#
# ### 3\. ECD 断裂状态模拟模块 (Scenarios)
#
# 这一部分是代码中最具生物学意义的地方。它模拟了：“ECD 能量打下去，到底打断了哪些二硫键？”
#
#   * **`get_cleavage_scenarios(pattern)`**:
#       * **背景**: 假设当前连接方式有 3 对二硫键 (A, B, C)。
#       * **逻辑**: 生成“幂集” (Power Set)。
#       * 它会列出所有可能的情况：
#           * 全断：[] (线性)
#           * 留1对：[A], [B], [C]
#           * 留2对：[A,B], [A,C], [B,C]
#           * 全没断：[A,B,C]
#       * 这覆盖了实验中可能发生的所有“部分还原”状态。
#
# -----
#
# ### 4\. 核心判定与计算模块 (判断 c/z 离子是否可见)
#
# 这是代码的“大脑”，决定了 Excel 表格中哪些是数值，哪些是 "SILENT"。
#
#   * **`is_bridged(cleavage_site, intact_bonds)`**:
#
#       * **作用**: 判断“沉默区域”。
#       * **原理**: 假设骨架在第 7 位断开 ($c_7 / z_{25}$)。如果此时存在一个**未断裂**的二硫键连接 C4 和 C10。
#       * **判定**: 因为 $4 < 7 < 10$，骨架虽然断了，但两头被 C4-C10 锁死，碎片飞不开。
#       * **输出**: 返回 `True`，后续代码会将该离子标记为 **SILENT**（不可见）。
#
#   * **`calculate_fragment_mass(...)`**:
#
#       * **作用**: 计算精确质量。
#       * **修正逻辑**:
#         1.  先算基础序列质量（氨基酸加和）。
#         2.  加上端基（c离子加 $NH_3$，z离子加 $OH-NH$）。
#         3.  **动态修正**: 检查这个碎片内部包含了几个完整的二硫键 (`internal_ss_count`)。每有一个，质量减去 `2.0156 Da`。
#
# -----
#
# ### 5\. 主循环与数据导出 (Main Loop)
#
# 函数 `analyze_sequence` 将上述步骤串联起来：
#
# 1.  **第一层循环**: 遍历 15 种连接模式 (Pattern 1\~15)。
# 2.  **第二层循环**: 遍历每种模式下的断裂状态 (全断、部分断、全不断)。
# 3.  **第三层循环**: 沿着肽链从头走到尾 (Position 1\~31)。
#       * 计算 **c 离子**：调用 `is_bridged` 判断是否可见 -\> 计算质量。
#       * 计算 **z 离子**：同理。
#       * 记录数据：将结果存入字典列表。
# 4.  **Pandas 导出**: 最后使用 `pd.DataFrame(data).to_excel` 将几千行数据瞬间写入 Excel。
#
# -----
#
# ### 代码逻辑流程图解
#
# ```text
# 输入序列 (SEQUENCE)
#        ↓
# 找到所有 Cys 位置 (例如 6个)
#        ↓
# 生成 15 种二硫键连接方式 (Patterns)
#        ↓
# [循环 1] 对于每一种连接方式 (例如: C4-C10, C14-C18, C23-C29):
#        ↓
#     生成断裂场景 (Scenarios): 比如 "C4-C10 没断，其他都断了"
#        ↓
#     [循环 2] 对于每一种场景:
#            ↓
#         [循环 3] 遍历每个氨基酸位置 (i = 1 to 31):
#                ↓
#             1. 模拟骨架断裂
#             2. 检查 is_bridged?
#                -> 是: 标记 "SILENT" (因为被拉住了)
#                -> 否: 计算质量 Mass = 序列重 - (内部完整二硫键数量 * 2H)
#             3. 存入列表
#        ↓
# 导出 Excel
# ```
#
# ### 总结
#
# 这段代码不仅仅是一个计算器，它是一个**虚拟实验平台**。它帮你穷举了所有理论上可能发生的化学事件，让你拿着生成的 Excel 表格，去和真实的质谱数据（Mass List）进行比对：
#
#   * 如果 Excel 说某段是 "SILENT"，而你没看到峰，说明模型匹配。
#   * 如果 Excel 说质量应该减 2 Da，而你的峰确实偏了 2 Da，说明该处有一个二硫键未断。

import itertools
import pandas as pd
import os

# ================= 用户配置区域 =================
# 目标多肽序列
SEQUENCE = "GCNILQPYWGCGRDFECLEECLMDSQYYQ"
OUTPUT_FILENAME = "ECD_Fragmentation_Analysis.xlsx"
# ==============================================

# --- 1. 基础质量常数 (单一同位素) ---
AA_MASS = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694, 'C': 103.00919,
    'E': 129.04259, 'Q': 128.05858, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
    'L': 113.08406, 'K': 128.09496, 'M': 131.04049, 'F': 147.06841, 'P': 97.05276,
    'S': 87.03203, 'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}

CONSTANTS = {
    'H_ATOM': 1.00783,
    'PROTON': 1.00728,
    'NH3': 17.02655,
    'OH': 17.00274,
    'NH': 15.0109,
    'SS_LOSS': 2 * 1.00783  # 形成一个二硫键失去的质量 (2H)
}


# --- 2. 核心计算逻辑函数 ---

def get_cys_indices(seq):
    """返回序列中所有 Cys 的位置 (1-based)"""
    return [i + 1 for i, char in enumerate(seq) if char == 'C']


def generate_disulfide_patterns(cys_list):
    """递归生成所有可能的二硫键配对方式"""
    if not cys_list:
        yield []
        return

    first = cys_list[0]
    rest = cys_list[1:]

    for i, partner in enumerate(rest):
        pair = (first, partner)
        remaining = rest[:i] + rest[i + 1:]
        for sub_pattern in generate_disulfide_patterns(remaining):
            yield [pair] + sub_pattern


def get_cleavage_scenarios(pattern):
    """
    对于一个给定的连接方式（如3对二硫键），
    生成所有可能的断裂状态（Power Set）。
    返回列表：每个元素是一个 list，包含当前【保持完整】的二硫键。
    """
    scenarios = []
    # pattern 的长度 (例如 3)
    n = len(pattern)
    # 生成 0 到 n 的所有组合长度
    for r in range(n + 1):
        # 找出所有保持完整的组合
        for intact_subset in itertools.combinations(pattern, r):
            scenarios.append(list(intact_subset))
    return scenarios


def is_bridged(cleavage_site, intact_bonds):
    """
    判断骨架断裂点是否位于某个完整的二硫键环内。
    cleavage_site: 氨基酸位置 i (表示断裂在 i 和 i+1 之间)
    """
    for c1, c2 in intact_bonds:
        start, end = sorted((c1, c2))
        # 如果断裂点在 start 和 end 之间
        if start <= cleavage_site < end:
            return True, f"C{start}-C{end}"
    return False, None


def calculate_fragment_mass(seq_fragment, internal_ss_count, type='c'):
    """计算片段质量，考虑内部二硫键引起的质量亏损"""
    base_mass = sum(AA_MASS[aa] for aa in seq_fragment)

    if type == 'c':
        # c-ion: Seq + NH3 + H+
        mass = base_mass + CONSTANTS['NH3'] + CONSTANTS['PROTON'] + CONSTANTS['H_ATOM']  # N-term H
    else:
        # z-dot ion: Seq + OH - NH + H + H+
        mass = base_mass + CONSTANTS['OH'] - CONSTANTS['NH'] + CONSTANTS['H_ATOM'] + CONSTANTS['PROTON']

    # 减去保留的二硫键质量 (每个键少 2H)
    mass -= internal_ss_count * CONSTANTS['SS_LOSS']
    return mass


def analyze_sequence(seq):
    cys_indices = get_cys_indices(seq)
    if len(cys_indices) % 2 != 0:
        print("警告：Cys数量为奇数，本程序仅处理完全配对情况。")
        return []

    all_patterns = list(generate_disulfide_patterns(cys_indices))
    results = []

    print(f"序列长度: {len(seq)}")
    print(f"Cys 位置: {cys_indices}")
    print(f"发现 {len(all_patterns)} 种连接模式。开始计算...")

    # 遍历每一种连接模式 (Pattern 1 ~ 15)
    for pat_idx, pattern in enumerate(all_patterns):
        pattern_str = "; ".join([f"C{a}-C{b}" for a, b in pattern])

        # 遍历该模式下，所有可能的 ECD 断裂程度 (全断、断1个、断2个...)
        # scenarios 是【保持完整】的二硫键列表
        scenarios = get_cleavage_scenarios(pattern)

        for scenario_idx, intact_bonds in enumerate(scenarios):
            # 描述当前状态
            if not intact_bonds:
                state_desc = "Fully Reduced (Linear)"
            elif len(intact_bonds) == len(pattern):
                state_desc = "Fully Oxidized (No S-S Broken)"
            else:
                remaining_str = ", ".join([f"C{a}-C{b}" for a, b in intact_bonds])
                state_desc = f"Partially Reduced (Intact: {remaining_str})"

            # --- 计算 c 离子 ---
            for i in range(len(seq) - 1):
                pos = i + 1  # 1-based
                frag_seq = seq[:pos]

                # 1. 检查是否沉默 (Silent)
                bridged, bridge_info = is_bridged(pos, intact_bonds)

                # 2. 计算内部有多少个完整的二硫键 (影响质量)
                internal_ss = 0
                for c1, c2 in intact_bonds:
                    if 1 <= c1 <= pos and 1 <= c2 <= pos:
                        internal_ss += 1

                if bridged:
                    mass = 0
                    status = "SILENT"
                    note = f"Bridged by {bridge_info}"
                else:
                    mass = calculate_fragment_mass(frag_seq, internal_ss, 'c')
                    status = "Observed"
                    note = f"Contains {internal_ss} S-S bonds" if internal_ss > 0 else ""

                results.append({
                    "Pattern_ID": pat_idx + 1,
                    "Connectivity": pattern_str,
                    "ECD_State": state_desc,
                    "Ion_Type": "c",
                    "Position": pos,
                    "Fragment": f"c{pos}",
                    "Sequence": frag_seq,
                    "Theo_Mass_MH+": mass if not bridged else "N/A",
                    "Status": status,
                    "Note": note
                })

            # --- 计算 z 离子 ---
            for i in range(len(seq) - 1):
                pos = i + 1  # 断裂在 pos 后
                z_pos = len(seq) - pos
                frag_seq = seq[pos:]

                # 1. 检查是否沉默
                bridged, bridge_info = is_bridged(pos, intact_bonds)

                # 2. 计算内部完整二硫键
                internal_ss = 0
                for c1, c2 in intact_bonds:
                    # z片段范围是 pos+1 到 len(seq)
                    if (pos + 1) <= c1 <= len(seq) and (pos + 1) <= c2 <= len(seq):
                        internal_ss += 1

                if bridged:
                    mass = 0
                    status = "SILENT"
                    note = f"Bridged by {bridge_info}"
                else:
                    mass = calculate_fragment_mass(frag_seq, internal_ss, 'z')
                    status = "Observed"
                    note = f"Contains {internal_ss} S-S bonds" if internal_ss > 0 else ""

                results.append({
                    "Pattern_ID": pat_idx + 1,
                    "Connectivity": pattern_str,
                    "ECD_State": state_desc,
                    "Ion_Type": "z·",  # z-dot
                    "Position": z_pos,
                    "Fragment": f"z{z_pos}",
                    "Sequence": frag_seq,
                    "Theo_Mass_MH+": mass if not bridged else "N/A",
                    "Status": status,
                    "Note": note
                })

    return results


# ================= 主程序执行 =================

if __name__ == "__main__":
    data = analyze_sequence(SEQUENCE)

    if data:
        print(f"计算完成，共生成 {len(data)} 条数据。正在写入 Excel...")
        df = pd.DataFrame(data)

        # 调整列顺序
        cols = ["Pattern_ID", "Connectivity", "ECD_State", "Ion_Type", "Fragment",
                "Theo_Mass_MH+", "Status", "Note", "Sequence"]
        df = df[cols]

        try:
            df.to_excel(OUTPUT_FILENAME, index=False)
            print(f"成功！文件已保存为: {os.path.abspath(OUTPUT_FILENAME)}")
            print("\nExcel 列说明:")
            print("1. Connectivity: 假设的初始二硫键连接方式 (共15种)")
            print("2. ECD_State: ECD打碎了哪些键，保留了哪些键 (Intact)")
            print("3. Status: 'Observed' 表示可见，'SILENT' 表示被二硫键拉住不可见")
            print("4. Note: 解释为什么不可见，或指出质量发生了多少个二硫键的位移 (-2Da per bond)")
        except Exception as e:
            print(f"写入Excel失败: {e}")
            print("请确保文件未被其他程序占用。")

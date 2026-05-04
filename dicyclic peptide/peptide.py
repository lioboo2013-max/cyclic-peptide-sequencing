import itertools
import pandas as pd
import os

# ================= 用户配置区域 =================
# 目标多肽序列
SEQUENCE = "AHTCPKLICYHHCFRECDSWQACGGGTTCLIP"
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
#!/usr/bin/env python3
"""
seq_alignment.py

基于结构对齐结果提取蛋白质序列对齐信息：
1. 加载之前保存的 PyMOL session（alignment.pse）。
2. 针对每个目标结构（除参考结构 "ref" 外），采用距离阈值方法获取对齐区域：
   - 利用 “ref and name CA within 2.0 of <target>” 及其对应目标选择，
     得到在对齐后空间距离小于指定阈值（本例中为2.0 Å）的CA原子集合，
     以近似判断两者对齐的残基。
3. 对每个目标结构，从这两个选择中依次提取残基（使用 CA 原子的残基信息），
   并利用预定义的三字母代码到一字母代码的映射构造局部对齐序列。
4. 将对齐序列以 FASTA 格式输出到文件（文件名格式为 {target}_seq_align.fasta）。

输出文件：
    - 每个目标结构生成一个 FASTA 文件（如：target1_seq_align.fasta）
日志记录：
    - 详细过程记录在 seq_alignment.log 中。

注意：
    • 此方法使用距离判断获得对齐残基集合，适合在结构对齐（如 cealign）后进行
      近似残基对应关系的提取，但不能完全替代专门的对齐映射算法。
    • 如遇残基数目不匹配情况，以较少者为准进行配对，并在日志中提醒用户。
"""

import os
import logging

try:
    import pymol2
except ImportError:
    raise ImportError("pymol2 模块未安装，请检查 PyMOL 是否正确配置。")

# 配置日志记录
logging.basicConfig(
    filename='seq_alignment.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# 三字母代码转换到一字母代码的字典
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

def extract_aligned_sequence(cmd, target):
    """
    对于指定目标结构，通过距离阈值方法提取对齐区域，构建参考结构(ref)
    与目标结构(target)的局部序列对齐。
    
    参数：
        cmd: PyMOL 命令接口
        target: 目标结构对象名称（字符串）
    
    返回：
        (ref_seq_aligned, target_seq_aligned)：对应的对齐序列字符串；
        如果提取失败，则返回 (None, None)。
    """
    # 定义距离阈值（单位 Ångström），可根据实际情况调节
    distance_threshold = 2.0

    # 基于空间距离建立选择：
    # 选择在 ref 中距离目标结构 target 的 CA 原子小于阈值的残基
    ref_sel = f"ref and name CA within {distance_threshold} of {target}"
    # 选择在目标结构 target 中距离 ref 的 CA 原子小于阈值的残基
    target_sel = f"{target} and name CA within {distance_threshold} of ref"
    
    logging.info(f"为目标结构 {target} 创建选择: '{ref_sel}' 与 '{target_sel}'")
    
    try:
        cmd.select("ref_aligned", ref_sel)
        cmd.select("target_aligned", target_sel)
    except Exception as e:
        logging.error(f"创建选择失败: {e}")
        return None, None

    try:
        ref_model = cmd.get_model("ref_aligned")
        target_model = cmd.get_model("target_aligned")
    except Exception as e:
        logging.error(f"获取模型数据失败: {e}")
        return None, None
    
    if not ref_model.atom or not target_model.atom:
        logging.error("未能获取到对齐区域中 CA 原子的模型数据。")
        return None, None

    # 定义排序函数：按链和残基编号排序（注意residue编号通常为字符串，此处尝试转为整数排序）
    def sort_key(atom):
        try:
            return (atom.chain, int(atom.resi))
        except ValueError:
            return (atom.chain, atom.resi)
    
    ref_atoms = sorted(ref_model.atom, key=sort_key)
    target_atoms = sorted(target_model.atom, key=sort_key)

    n_aligned = min(len(ref_atoms), len(target_atoms))
    if n_aligned == 0:
        logging.error(f"目标结构 {target} 的对齐区域为空。")
        return None, None

    if len(ref_atoms) != len(target_atoms):
        logging.warning(
            f"目标结构 {target} 中对齐 CA 原子数目不匹配：ref={len(ref_atoms)} vs {target}={len(target_atoms)}。"
            " 取较小者进行配对。"
        )

    ref_seq_aligned = ""
    target_seq_aligned = ""
    for i in range(n_aligned):
        ref_resn = ref_atoms[i].resn.upper()
        target_resn = target_atoms[i].resn.upper()
        ref_letter = THREE_TO_ONE.get(ref_resn, "X")
        target_letter = THREE_TO_ONE.get(target_resn, "X")
        ref_seq_aligned += ref_letter
        target_seq_aligned += target_letter

    logging.info(f"{target}: 提取对齐残基数 = {n_aligned}")
    return ref_seq_aligned, target_seq_aligned

def main():
    logging.info("==== 开始提取蛋白质序列对齐信息流程 ====")
    session_file = "alignment.pse"

    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd

        # 加载已有的 PyMOL session
        logging.info(f"加载 PyMOL session 文件: {session_file}")
        try:
            cmd.load(session_file, "session")
        except Exception as e:
            logging.error(f"加载 session 失败: {e}")
            return

        # 获取所有对象名称，排除 "ref"（参考结构）和 session 标记
        objects = cmd.get_object_list()
        target_objects = [obj for obj in objects if obj not in ("ref", "session")]
        if not target_objects:
            logging.error("在 session 中未找到目标结构对象。")
            return

        for target in target_objects:
            logging.info(f"开始处理目标结构: {target}")
            ref_aln, tgt_aln = extract_aligned_sequence(cmd, target)
            if ref_aln is None or tgt_aln is None:
                logging.error(f"目标结构 {target} 的对齐序列提取失败。")
                continue

            # 生成 FASTA 格式的对齐结果
            fasta_output = (
                f">ref_aligned\n{ref_aln}\n"
                f">{target}_aligned\n{tgt_aln}\n"
            )
            output_filename = f"{target}_seq_align.fasta"
            try:
                with open(output_filename, "w") as outfile:
                    outfile.write(fasta_output)
                logging.info(f"写入对齐序列文件: {output_filename}")
            except Exception as e:
                logging.error(f"写入文件 {output_filename} 失败: {e}")

    logging.info("==== 蛋白质序列对齐信息提取流程结束 ====")

if __name__ == "__main__":
    main()

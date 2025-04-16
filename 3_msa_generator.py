#!/usr/bin/env python3
"""
msa_generator.py

该脚本用于将脚本二生成的各个局部蛋白序列对齐文件（FASTA格式）合并为一份全局的多序列比对（MSA）文件。
假定每个局部对齐文件包含两条序列：
    (1) 参考序列，标题为 "ref_aligned"
    (2) 目标结构的对齐序列，标题格式为 "{target}_aligned"

所有局部对齐均来源于参考结构相同区域，因此参考序列应保持一致（若不一致，则发出警告）。
最终将合并结果输出为 aligned.msa.fasta 文件，以供后续分析使用。
"""

import os
import glob
import logging

# 配置日志记录
logging.basicConfig(
    filename='msa_generator.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def parse_fasta(filename):
    """
    解析 FASTA 文件，返回字典，键为标题（去除 '>'），值为全部拼接的序列字符串。
    """
    sequences = {}
    with open(filename, 'r') as f:
        lines = f.read().splitlines()
    
    header = None
    seq_lines = []
    for line in lines:
        if line.startswith(">"):
            if header is not None:
                sequences[header] = ''.join(seq_lines)
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        sequences[header] = ''.join(seq_lines)
    return sequences

def main():
    logging.info("==== 开始合并序列对齐生成 MSA 文件流程 ====")
    
    # 搜索所有局部对齐 FASTA 文件，文件命名格式必须为 *_seq_align.fasta
    fasta_files = glob.glob("*_seq_align.fasta")
    if not fasta_files:
        logging.error("未找到任何局部对齐 FASTA 文件。")
        return
    
    # 保存所有参考序列与目标序列
    ref_sequences = []  # 用于存储每个文件中提取到的参考序列
    targets = {}        # 保存目标结构对齐序列，键为目标结构名称（不含后缀 _aligned）
    
    # 遍历所有文件进行解析
    for fasta_file in fasta_files:
        logging.info(f"处理文件: {fasta_file}")
        try:
            seq_dict = parse_fasta(fasta_file)
        except Exception as e:
            logging.error(f"解析文件 {fasta_file} 失败: {e}")
            continue
        
        # 检查文件中必须包含 reference 对齐序列（标题为 "ref_aligned"）以及目标结构对齐序列
        if "ref_aligned" not in seq_dict:
            logging.error(f"文件 {fasta_file} 缺少 'ref_aligned' 序列。")
            continue
        
        # 找出除 "ref_aligned" 之外的唯一目标序列
        target_keys = [key for key in seq_dict if key != "ref_aligned"]
        if len(target_keys) != 1:
            logging.error(f"文件 {fasta_file} 格式不符合预期，目标序列数量异常：{target_keys}")
            continue
        
        ref_seq = seq_dict["ref_aligned"]
        target_key = target_keys[0]
        target_seq = seq_dict[target_key]
        
        ref_sequences.append(ref_seq)
        # 去除后缀 "_aligned" 得到目标的名称
        target_name = target_key.replace("_aligned", "")
        targets[target_name] = target_seq

    if not ref_sequences:
        logging.error("没有有效的参考序列数据，终止处理。")
        return

    # 检查所有参考序列是否一致
    first_ref = ref_sequences[0]
    consistent = all(seq == first_ref for seq in ref_sequences)
    if not consistent:
        logging.warning("检测到不同的参考序列，可能局部对齐片段对应不同区域。")
    else:
        logging.info("所有参考序列一致。")
    
    # 构造最终 MSA 字典：参考序列及所有目标序列
    msa = {}
    msa["ref"] = first_ref
    for target_name, sequence in targets.items():
        msa[target_name] = sequence

    # 写入最终 MSA 文件，格式为 FASTA
    output_filename = "aligned.msa.fasta"
    try:
        with open(output_filename, "w") as outfile:
            for seq_name, sequence in msa.items():
                outfile.write(f">{seq_name}\n")
                outfile.write(f"{sequence}\n")
        logging.info(f"最终 MSA 写入文件：{output_filename}")
    except Exception as e:
        logging.error(f"写入 MSA 文件 {output_filename} 失败: {e}")

    logging.info("==== MSA 生成流程结束 ====")

if __name__ == "__main__":
    main()

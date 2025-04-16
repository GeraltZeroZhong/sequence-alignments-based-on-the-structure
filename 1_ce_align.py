#!/usr/bin/env python3
"""
ce_align.py

该脚本用于对 /pdbs 文件夹内的所有 .pdb 文件（除 ref.pdb 外）
进行结构批量对齐。以 ref.pdb 为参考结构，利用 PyMOL 内置的
Combinatorial Extension (CE) 对齐算法（cealign 命令）进行比对，
记录每次比对的 RMSD 值，并保存对齐日志以及 PyMOL session 文件。

输出文件：
    - alignment.log       日志文件：记录比对过程及各结构的 RMSD 值
    - alignment.pse       保存包含所有结构及对齐结果的 PyMOL session
"""

import os
import glob
import logging
from datetime import datetime

import os

# 采用 pymol2 模块以编程方式启动 PyMOL 环境，注意确保安装了 PyMOL2
try:
    import pymol2
except ImportError:
    raise ImportError("未能导入 pymol2 模块，请检查 PyMOL 是否正确安装。")

# 配置日志输出
logging.basicConfig(
    filename='alignment.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def main():
    logging.info("==== 开始 CE 对齐流程 ====")
    
    # 设置 pdb 文件所在目录和参考结构文件
    pdb_folder = os.path.join(os.getcwd(), "pdbs")
    ref_filename = "ref.pdb"
    ref_filepath = os.path.join(pdb_folder, ref_filename)
    
    # 检查参考结构是否存在
    if not os.path.exists(ref_filepath):
        logging.error(f"参考结构 {ref_filepath} 不存在！请检查文件路径。")
        return
    
    # 查找目录中的所有 .pdb 文件，并剔除参考结构文件
    pdb_files = glob.glob(os.path.join(pdb_folder, "*.pdb"))
    pdb_files = [f for f in pdb_files if os.path.basename(f) != ref_filename]
    
    if not pdb_files:
        logging.error(f"在目录 {pdb_folder} 中未找到除 {ref_filename} 外的其他 .pdb 文件！")
        return

    # 使用 pymol2 启动 PyMOL 环境
    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd  # 获取 PyMOL 的命令接口
        
        # 加载参考结构
        logging.info(f"加载参考结构：{ref_filepath}")
        cmd.load(ref_filepath, "ref")
        
        # 遍历其他 pdb 文件，依次加载并比对
        for pdb_path in pdb_files:
            base_name = os.path.splitext(os.path.basename(pdb_path))[0]
            logging.info(f"加载目标结构：{pdb_path}")
            cmd.load(pdb_path, base_name)
            
            # 执行 CE 对齐，比对参考结构 "ref" 与目标结构 base_name
            logging.info(f"进行 cealign 对齐：参考结构 'ref' 与 '{base_name}'")
            try:
                result = cmd.cealign("ref", base_name)
            except Exception as e:
                logging.error(f"对齐 {base_name} 时出现异常: {e}")
                continue
            
            # 提取 RMSD 值，PyMOL 返回的 result 可能为字典或其他数据类型
            rmsd = None
            if isinstance(result, dict):
                # 可能的 key 为 "RMS" 或 "RMSD"
                rmsd = result.get("RMS", result.get("RMSD", None))
            else:
                # 当返回值非字典时，可尝试直接转换
                try:
                    rmsd = float(result)
                except Exception:
                    rmsd = "N/A"

            logging.info(f"对齐结果 {base_name}：RMSD = {rmsd}")
        
        # 保存 PyMOL session，包含所有加载及对齐后的结构
        session_file = "alignment.pse"
        logging.info(f"保存 PyMOL session 为 {session_file}")
        cmd.save(session_file)
        
    logging.info("==== CE 对齐流程结束 ====")

if __name__ == "__main__":
    main()
